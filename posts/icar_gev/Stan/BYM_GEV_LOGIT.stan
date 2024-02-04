
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
  real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
    return - 0.5 * dot_self((phi[node1] - phi[node2])) 
      + normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
}

data {
  int<lower = 0> N_years;
  int<lower = 0> N_stations;
  matrix[N_years, N_stations] precip;


  int<lower = 0> N_neighbors;
  array[N_neighbors] int node1;
  array[N_neighbors] int node2;
  real<lower = 0> scaling_factor;
}



parameters {
  vector[N_stations] psi_random;
  vector[N_stations] psi_spatial;
  real<lower = 0> sigma_psi;
  real mu_psi;
  real logit_rho_psi;
  
  vector[N_stations] tau_random;
  vector[N_stations] tau_spatial;
  real<lower = 0> sigma_tau;
  real mu_tau;
  real logit_rho_tau;
  
  vector[N_stations] phi_random;
  vector[N_stations] phi_spatial;
  real<lower = 0> sigma_phi;
  real mu_phi;
  real logit_rho_phi;
  
  vector[N_stations] gamma_random;
  vector[N_stations] gamma_spatial;
  real<lower = 0> sigma_gamma;
  real mu_gamma;
  real logit_rho_gamma;
}

transformed parameters {
  real<lower = 0, upper = 1> rho_psi = inv_logit(logit_rho_psi);
  real<lower = 0, upper = 1> rho_tau = inv_logit(logit_rho_tau);
  real<lower = 0, upper = 1> rho_phi = inv_logit(logit_rho_phi);
  real<lower = 0, upper = 1> rho_gamma = inv_logit(logit_rho_gamma);
  
  vector[N_stations] psi = mu_psi + sigma_psi * (sqrt(rho_psi / scaling_factor) * psi_spatial + sqrt(1 - rho_psi) * psi_random);
  vector[N_stations] tau = mu_tau + sigma_tau * (sqrt(rho_tau / scaling_factor) * tau_spatial + sqrt(1 - rho_tau) * tau_random);
  vector[N_stations] phi = mu_phi + sigma_phi * (sqrt(rho_phi / scaling_factor) * phi_spatial + sqrt(1 - rho_phi) * phi_random);
  vector[N_stations] gamma = mu_gamma + sigma_gamma * (sqrt(rho_gamma / scaling_factor) * gamma_spatial + sqrt(1 - rho_gamma) * gamma_random);
  
  vector<lower = 0>[N_stations] mu0 = exp(psi);
  vector<lower = 0>[N_stations] sigma = exp(psi + tau);
  vector<lower = -0.5, upper = 0.5>[N_stations] xi = inv_logit(phi) - 0.5;
  vector<lower = -0.01, upper = 0.01>[N_stations] delta = 0.02 * inv_logit(gamma) - 0.01;
}

model {
  for (i in 1:N_stations) {
    precip[ , i] ~ gevt(mu0[i], sigma[i], xi[i], delta[i]);
  }

  psi_spatial ~ icar_normal(N_neighbors, node1, node2);
  psi_random ~ std_normal();
  sigma_psi ~ exponential(1);
  logit_rho_psi ~ std_normal();
  mu_psi ~ normal(2.2, 1);
  
  tau_spatial ~ icar_normal(N_neighbors, node1, node2);
  tau_random ~ std_normal();
  sigma_tau ~ exponential(1);
  logit_rho_tau ~ std_normal();
  mu_tau ~ normal(-0.9, 1);
  
  phi_spatial ~ icar_normal(N_neighbors, node1, node2);
  phi_random ~ std_normal();
  sigma_phi ~ exponential(1);
  logit_rho_phi ~ std_normal();
  mu_phi ~ normal(0, 1);
  
  gamma_spatial ~ icar_normal(N_neighbors, node1, node2);
  gamma_random ~ std_normal();
  sigma_gamma ~ exponential(1);
  logit_rho_gamma ~ std_normal();
  mu_gamma ~ normal(0, 1);

}

