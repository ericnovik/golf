functions {
  real theta0(real x, real R, real r) {
    return asin((R - r) / x);
  }
}

data {
  int N;
  int<lower = 0> tries[N];
  int<lower = 0> successes[N];
  real<lower = 0> dist[N];
}

transformed data {
  real R;
  real r;
  R <- 4.25 / 2;
  r <- 1.68 / 2;
}

parameters {
  real<lower = 0> sigma;
}

model {
  real p[N];
  
  for (n in 1:N) {
    p[n] <- 2 * Phi(theta0(dist[n], R, r) / sigma) - 1;
  }
  
  sigma ~ cauchy(0, 2.5);
  successes ~ binomial(tries, p);
}

generated quantities {
  real sigma_degrees;
  sigma_degrees <- 180/pi() * sigma;
}
