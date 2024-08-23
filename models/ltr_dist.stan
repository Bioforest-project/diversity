data {
  int<lower=0> N ; // observation #
  int<lower=0> P ; // plot #
  vector[N] diversity ; // diversity values
  vector[N] time ; // time values
  vector[P] disturbance ; // disturbance values
  int<lower=1, upper=P> plot[N] ; // plot index
}
parameters {
  vector<upper=0>[P] a_p ; // plot quadratic parameter
  vector[P] b_p ; // plot linear parameter
  vector[P] c_p ; // plot intercept parameter
  real d ; // disturbance slope parameter
  real d_max0 ; // disturbance intercept parameter
  real<lower=0> sigma ; // diversity variation
  real<lower=0> sigma_d ; // disturbance variation
}
transformed parameters {
  vector[P] t_opt = - b_p ./ (2*a_p) ; // optimum time
  vector[P] d_max = a_p .* t_opt .* t_opt + b_p .* t_opt + c_p ; // maximum diversity
}
model {
  diversity ~ normal(a_p[plot] .* time .* time + b_p[plot] .* time + c_p[plot], sigma) ;
  d_max ~ normal(d_max0 + d*disturbance, sigma_d) ;
}
