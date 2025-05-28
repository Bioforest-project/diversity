data {
  int<lower=0> n_obs;
  int<lower=0> n_site;
  int<lower=0> n_plot;
  vector[n_obs] y;
  array[n_obs] int<lower=0, upper=n_site> site;
  array[n_obs] int<lower=0, upper=n_plot> plot;
}
parameters {
  vector<lower=0>[n_site] mu;
  vector[n_plot] delta;
  real<lower=0> sigma;
  real<lower=0> sigma_p;
}
model {
  y ~ normal(mu[site] + delta[plot], sigma);
  delta ~ normal(0, sigma_p);
  // delta ~ lognormal(0, sigma_p);
  // log(y) ~ normal(log(mu[site]) .* delta[plot], sigma);
  // delta ~ lognormal(0, sigma_p);
}
