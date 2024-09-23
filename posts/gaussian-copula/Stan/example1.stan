functions {
  real gaussian_copula_lpdf(vector u, matrix L) {
    int D = num_elements(u);
    vector[D] z = inv_Phi(u);
    return multi_normal_cholesky_lpdf(z | rep_vector(0, D), L) - normal_lpdf(z | 0, 1);
  }

  real exponential_icdf(real u, real lambda) {
    return -log(1 - u) / lambda;
  }

}
data {
  int<lower = 0> N;
  int<lower = 0> D;
  matrix[N, D] X;
}

parameters {
  vector[D] lambda;
  cholesky_factor_corr[D] L;
}

model {
  matrix[N, D] U;
  for (i in 1:N) {
    for (j in 1:D) {
      target += exponential_lpdf(X[i, j] | lambda[j]);
      U[i, j] = exponential_cdf(X[i, j] | lambda[j]);
    }
    target += gaussian_copula_lpdf(to_vector(U[i, ]) | L);
  }
  
  target += lkj_corr_cholesky_lpdf(L | 1.0);
}

generated quantities {
  corr_matrix[D] Sigma = multiply_lower_tri_self_transpose(L);
  matrix[N, D] yrep;

  {
    matrix[N, D] Z_rep;
    matrix[N, D] U_rep;

    for (i in 1:N) {
      Z_rep[i, ] = to_row_vector(multi_normal_cholesky_rng(rep_vector(0, D), L));
      for (j in 1:D) {
        U_rep[i, j] = Phi(Z_rep[i, j]);
        yrep[i, j] = exponential_icdf(U_rep[i, j], lambda[j]);
      }
    }
  }
}
