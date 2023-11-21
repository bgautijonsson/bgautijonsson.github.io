
functions{
  /*
  Fit a GEV distribution with trend
  */
  real gevt_lpdf(vector y, real mu0, real sigma, real xi, real delta) {
    int N = rows(y);
    vector[N] z;
    vector[N] mu;
    real lp = 0;

    for (i in 1:N) {
      mu[i] = mu0 * (1 + delta * (i - 1));


      if (z[i] <= -1) {
        reject("found incompatible variable values");
      }
      
      if (abs(xi) < 1e-12) {
        z[i] = (y[i] - mu[i]) / sigma;
        lp += - z[i] - exp(-z[i]);
      } else {
        z[i] = xi * (y[i] - mu[i]) / sigma;
        lp += - (1 + 1 / xi) * log(1 + z[i]) - pow((1 + z[i]), -1/xi);
      }
    }

    lp += -N * log(sigma);

    return lp;
  }
  /*
  Put an ICAR prior on coefficients
  */
  real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2, real sigma) {
    return - N * log(sigma) - 0.5 * dot_self((phi[node1] - phi[node2]) / sigma);
  }
}

data {
  int<lower = 0> N_years;
  int<lower = 0> N_stations;
  matrix[N_years, N_stations] precip;


  int<lower = 0> N_neighbors;
  array[N_neighbors] int node1;
  array[N_neighbors] int node2;
}



parameters {
  vector[N_stations] psi;
  real<lower = 0> sigma_psi;
  
  vector[N_stations] tau;
  real<lower = 0> sigma_tau;
  
  vector[N_stations] phi;
  real<lower = 0> sigma_phi;
  
  vector[N_stations] gamma;
  real<lower = 0> sigma_gamma;
  
}

transformed parameters {
  vector<lower = 0>[N_stations] mu0 = exp(psi);
  vector<lower = 0>[N_stations] sigma = exp(psi + tau);
  vector<lower = -0.5>[N_stations] xi = exp(phi) - 0.5;
  vector<lower = -0.008, upper = 0.008>[N_stations] delta = 0.008 * inv_logit(gamma);
}

model {
  for (i in 1:N_stations) {
    precip[ , i] ~ gevt(mu0[i], sigma[i], xi[i], delta[i]);
  }

  psi ~ icar_normal(N_neighbors, node1, node2, sigma_psi);
  tau ~ icar_normal(N_neighbors, node1, node2, sigma_tau);
  phi ~ icar_normal(N_neighbors, node1, node2, sigma_phi);
  gamma ~ icar_normal(N_neighbors, node1, node2, sigma_gamma);

}

