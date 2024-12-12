functions {

real student_t_icdf(real p, real nu) {
  real eps = machine_precision();
  real d_epsilon = machine_precision();
  real d_max = positive_infinity();
  real d_min = machine_precision();
  int d_mant_dig = 53;
  real q;
  
  if (is_nan(p) || is_nan(nu)) {
    return p + nu;
  }
  
  if (nu <= 0) {
    reject("Invalid value for nu");
  }
  
  if (nu < 1) {
    // Find the upper and lower bounds
    real accu = machine_precision();
    real Eps = machine_precision();

    if (p > 1 - d_epsilon) {
      return positive_infinity();
    }

    real pp = min({1 - d_epsilon, p * (1 + Eps)});
    real ux = 1;
    while (student_t_cdf(ux | nu, 0., 1.) < pp) {
      ux = ux * 2;
    }
    pp = p * (1 - Eps);
    real lx = -1;
    while (student_t_cdf(lx | nu, 0., 1.) > pp) {
      lx = lx * 2;
    }
    
    // Find the quantile using interval halving
    real nx = 0.5 * (lx + ux);
    int iter = 0;
    while ((ux - lx) / abs(nx) > accu && iter < 1000) {
      iter += 1;
      if (student_t_cdf(nx | nu, 0., 1.) > p) {
        ux = nx;
      } else {
        lx = nx;
      }
      nx = 0.5 * (lx + ux);
    }
    return 0.5 * (lx + ux);
  }
  
  if (nu > 1e20) {
    return inv_Phi(p);
  }
  
  int neg = p < 0.5 ? 1 : 0;
  int is_neg_lower = neg;
  real P = neg == 1 ? 2 * p : 2 * (0.5 - p + 0.5);
  
  P = min({max({P, 0}), 1});
  
  if (abs(nu - 2) < eps) {
    if (P > d_min) {
      if (3 * P < d_epsilon) {
        q = 1 / sqrt(P);
      } else if (P > 0.9) {
        q = (1 - P) * sqrt(2 / (P * (2 - P)));
      } else {
        q = sqrt(2 / (P * (2 - P)) - 2);
      }
    } else {
      q = positive_infinity();
    }
  } else if (nu < 1. + eps) {
    if (P == 1.) {
      q = 0;
    } else if (P > 0) {
      q = 1 / tan(pi() * p / 2);
    } else {
      q = negative_infinity();
    }
  } else {
    real x = 0;
    real y;
    real log_P2 = 0;
    real a = 1 / (nu - 0.5);
    real b = 48 / (a * a);
    real c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
    real d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * pi() / 2) * nu;
    
    y = pow((d * P), (2.0 / nu));
    int P_ok = y >= d_epsilon ? 1 : 0;
    
    if (P_ok != 1) {
      log_P2 = is_neg_lower == 1 ? log(p) : log1m_exp(p);
      x = (log(d) + log2() + log_P2) / nu;
      y = exp(2 * x);
    }
    
    if ((nu < 2.1 && P > 0.5) || y > 0.05 + a) {
      if (P_ok == 1) {
        x = inv_Phi(0.5 * P);
      } else {
        x = inv_Phi(log_P2);
      }
      
      y = square(x);
      
      if (nu < 5) {
        c += 0.3 * (nu - 4.5) * (x + 0.6);
      }
      
      c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
      y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1)
          * x;
      y = expm1(a * square(y));
      q = sqrt(nu * y);
    } else if (P_ok != 1 && x < -log2() * d_mant_dig) {
      q = sqrt(nu) * exp(-x);
    } else {
      y = ((1
            / (((nu + 6) / (nu * y) - 0.089 * d - 0.822) * (nu + 2) * 3)
            + 0.5 / (nu + 4))
           * y - 1)
          * (nu + 1) / (nu + 2) + 1 / y;
      
      q = sqrt(nu * y);
    }
    
    if (P_ok == 1) {
    int it = 0;
    while (it < 10) {
      y = exp(student_t_lpdf(q | nu, 0., 1.));
      if (y <= 0 || is_inf(y)) {
        break;
      }
      
      real t = (exp(student_t_lccdf(q | nu, 0., 1.)) - P / 2) / y;
      if (abs(t) <= 1e-14 * abs(q)) {
        break;
      }
      
      q = q + t * (1. + t * q * (nu + 1) / (2 * (q * q + nu)));
      
      it += 1;
    }
    }
  }
  
  if (neg) {
    q = -q;
  }
  
  return q;
}

vector student_t_icdf(vector u, real nu) {
  int N = num_elements(u);
  vector[N] x;

  for (i in 1 : N) {
    x[i] = student_t_icdf(u[i], nu);
  }
  
  return x;
}


  real t_copula_lpdf(vector u, matrix L, real nu) {
    int D = num_elements(u);
    vector[D] x;
    real logp;

    // Transform U to X via the inverse t CDF
    for (d in 1:D) {
      x[d] = student_t_icdf(u[d], nu);
    }

    // Multivariate t density
    // Stan provides: multi_student_t_cholesky_lpdf(x | nu, mu, L)
    // with mu a vector and L a cholesky factor of scale matrix.
    logp = multi_student_t_cholesky_lpdf(x | nu, rep_vector(0, D), L);

    // Subtract the product of marginal t densities:
    // Each marginal: student_t_lpdf(x[d]|nu,0,1)
    for (d in 1:D) {
      logp -= student_t_lpdf(x[d] | nu, 0, 1);
    }

    return logp;
  }

  real exponential_icdf(real u, real lambda) {
    return -log1m(u) / lambda;
  }
}

data {
  int<lower = 0> N;
  int<lower = 0> D;
  matrix[N, D] Y;
}

parameters {
  vector<lower=0>[D] lambda;
  cholesky_factor_corr[D] L;
  real<lower=0> nu;   // degrees of freedom for the t-copula
}

model {
  matrix[N, D] U;
  for (i in 1:N) {
    // Transform data to uniforms using exponential CDF
    for (j in 1:D) {
      target += exponential_lpdf(Y[i, j] | lambda[j]);
      U[i, j] = exponential_cdf(Y[i, j] | lambda[j]);
    }
    // Add the t-copula contribution
    target += t_copula_lpdf(to_vector(U[i, ]) | L, nu);
  }

  // Priors
  target += lkj_corr_cholesky_lpdf(L | 1.0);
  // Prior on nu can be chosen, e.g.:
  target += exponential_lpdf(nu | 1);
  target += exponential_lpdf(lambda | 1);
}

generated quantities {
  corr_matrix[D] Sigma = multiply_lower_tri_self_transpose(L);
  matrix[N, D] yrep;

  {
    matrix[N, D] U_rep;
    matrix[N, D] Z_rep;

    // Generate new samples:
    // 1. Draw Z_rep from multivariate t
    for (i in 1:N) {
      Z_rep[i] = (multi_student_t_cholesky_rng(nu, rep_vector(0, D), L))';
      // Convert from t-variate to uniforms using t CDF
      for (j in 1:D) {
        U_rep[i, j] = student_t_cdf(Z_rep[i, j] | nu, 0, 1);
        // Convert uniforms to exponential
        yrep[i, j] = exponential_icdf(U_rep[i, j], lambda[j]);
      }
    }
  }
}
