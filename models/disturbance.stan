data {
  int<lower=0> n; // # obs
  int<lower=0> s; // # sites
  int<lower=0> p; // # plots
  vector[n] ba; // basal area
  vector[n] y; // forest diversity
  vector[n] t; // time
  array[n] int<lower=0, upper=n> site; // site index
  array[n] int<lower=0, upper=n> plot; // plot index
  vector[p] dist_p; // disturbance
  array[p] int<lower=0, upper=n> site_plot; // site index in plots
  vector[s] mu_theta_s; // mean site asymptotic ba in controls
  vector[s] sigma_theta_s; // sd site asymptotic ba in controls
}
transformed data {
  vector[55] t_pred;
  for(i in 1:55)
    t_pred[i] = i;
}
parameters {
  vector<lower=0.1, upper=1>[p] phi_p;
  vector<lower=0, upper=0.5>[p] lambda_p;
  real<lower=0, upper=0.5> mu_lambda;
  real<lower=0> sigma_lambda;
  vector<lower=5-3, upper=40-3>[s] tau0_s;
  real<lower=5-3, upper=40-3> mu_tau0;
  real<lower=0> sigma_tau;
  vector<lower=(mu_theta_s-sigma_theta_s)*0.1, upper=(mu_theta_s+sigma_theta_s)*10>[s] theta_s;
  real<lower=0> sigma;
  vector<lower=-phi_p, upper=2>[p] delta_p;
  real<lower=0> sigma_delta;
  vector<lower=-4, upper=4>[s] delta0_s;
  vector<lower=-4, upper=4>[s] gamma_s;
}
transformed parameters {
  vector[n] ltp = 1 - exp(-lambda_p[plot] .* t);
  vector[n] stp = delta_p[plot] .*
                  (t ./ tau0_s[site] .* exp(1 - t ./ tau0_s[site])) .*
                  (t ./ tau0_s[site] .* exp(1 - t ./ tau0_s[site]));
  vector[n] mu = theta_s[site] .* (phi_p[plot] + ltp.*(1 - phi_p[plot]) + stp);
}
model {
  log(y) ~ normal(log(mu), sigma);
  lambda_p ~ normal(mu_lambda, sigma_lambda);
  delta_p ~ cauchy(delta0_s[site_plot] + gamma_s[site_plot] .* dist_p, sigma_delta);
  tau0_s ~ normal(mu_tau0, sigma_tau);
  for(i in 1:s)
    theta_s[i] ~ normal(mu_theta_s[i], sigma_theta_s[i]);
  sigma ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
}
generated quantities {
  matrix[55, p] y_pred;
  for(i in 1:p)
    y_pred[,i] = theta_s[site_plot[i]] .* 
                  (phi_p[i] + 
                  (1 - exp(-lambda_p[i] * t_pred)) 
                  .* (1 - phi_p[i]) + 
                  (delta_p[i] * 
                 (t_pred / tau0_s[site_plot[i]] .* 
                 exp( 1 - t_pred ./ tau0_s[site_plot[i]] )) .*
                 (t_pred / tau0_s[site_plot[i]] .* 
                 exp( 1 - t_pred ./ tau0_s[site_plot[i]] ))));
}
