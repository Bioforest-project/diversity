data {
  int<lower=0> N ; // observation #
  int<lower=0> P ; // plot #
  vector[N] ba ; // basal area values
  vector[N] time ; // time values
  array[N] int<lower=1, upper=P> plot ; // plot index
}
parameters {
  vector<lower=0>[P] ba_post ; //  plot post disturbance basal area
  vector<lower=0>[P] delta_ba ; // plot disturbance intensity
  vector<lower=0>[P] alpha ; // plot decay parameter
  real<lower=0> sigma ; // basal area variation variation
}
model {
  ba ~ normal(ba_post[plot] + delta_ba[plot] .* exp(-alpha[plot] .* time), sigma) ;
  ba_post ~ normal(20, 5) ;
  delta_ba ~ normal(0, 5) ;
  alpha ~ normal(0, 1) ;
}
generated quantities {
  vector[P] ba_pre = ba_post + delta_ba ; // plot pre disturbance basal area
  vector[P] t_95 = 3/alpha ; // maximum time
}
