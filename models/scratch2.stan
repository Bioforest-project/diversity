data {
  int<lower=0> n_rec ; // obs reco
  int<lower=0> n_old ; // obs old
  int<lower=0> n_pre ; // obs prelog
  int<lower=0> n_site ;
  int<lower=0> n_plot ;
  vector[n_rec] stem_rec ;
  vector[n_old] stem_old ;
  vector[n_pre] stem_pre ;
  vector[n_rec] time ;
  int<lower=0, upper=n_site> site_rec[n_rec] ; 
  int<lower=0, upper=n_site> site_old[n_old] ; 
  int<lower=0, upper=n_site> site_pre[n_pre] ; 
  int<lower=0, upper=n_plot> plot_rec[n_rec] ; 
  int<lower=0, upper=n_plot> plot_old[n_old] ; 
  int<lower=0, upper=n_plot> plot_pre[n_pre] ; 
}
parameters {
  real<lower=0, upper=2000> mu_theta0 ; // starting point
  real<lower=0> sigma_theta0 ;
  vector<lower=0, upper=2000>[n_plot] theta0_p ;
  real<lower=0, upper=0.5> mu_lambda ; // recovery rate
  real<lower=0> sigma_lambda ;
  vector<lower=0, upper=0.5>[n_plot] lambda_p ;
  real<lower=0, upper=2000> mu_thetaInf ; // ending point
  real<lower=0> sigma_thetaInf ;
  vector<lower=0, upper=2000>[n_site] thetaInf_s ;
  real<lower=0> sigma_old ;
  real<lower=0> sigma_pre ;
  real<lower=0> sigma_rec ;
}
model {
  real mu_old[n_old] ;
  real mu_pre[n_pre] ;
  real mu_rec[n_rec] ;
  for(i in 1:n_old){
    mu_old[i] = thetaInf_s[site_old[i]] ;
  }
  for(j in 1:n_pre){
    mu_pre[j] = thetaInf_s[site_pre[j]] ;
  }
  for(k in 1:n_rec){
    mu_rec[k] = theta0_p[plot_rec[k]] + 
                (thetaInf_s[site_rec[k]] - theta0_p[plot_rec[k]]) *
                (1 - exp(-lambda_p[plot_rec[k]] * time[k])) ;
  }
  stem_old ~ lognormal(log(mu_old), sigma_old) ;
  stem_pre ~ lognormal(log(mu_pre), sigma_pre) ;
  stem_rec ~ lognormal(log(mu_rec), sigma_rec) ;
  theta0_p ~ lognormal(log(mu_theta0), sigma_theta0) ;
  thetaInf_s ~ lognormal(log(mu_thetaInf), sigma_thetaInf) ;
  lambda_p ~ lognormal(log(mu_lambda), sigma_lambda) ;
}
generated quantities {
  vector[n_rec] pred_stem_rec ;
  for(l in 1:n_rec){
    pred_stem_rec[l] = lognormal_rng(theta0_p[plot_rec[l]] * (1 - exp(-lambda_p[plot_rec[l]] * time[l])), sigma_rec) ;
  }
}
