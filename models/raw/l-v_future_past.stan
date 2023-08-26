functions {
    real[] dz_dt(real t, real[] z, real[] theta,
              real[] x_r, int[] x_i) {
              real u = z[1];
              real v = z[2];

              real alpha = theta[1];
              real beta = theta[2];
              real gamma = theta[3];
              real delta = theta[4];

              real du_dt = (alpha - beta * v) * u;
              real dv_dt = (-gamma + delta * u) * v;
              return { du_dt, dv_dt };
    }
}

data {
  int<lower = 0> N;          // number of observed measurements
  array[N] real ts;                // observed measurement times
  array[50] real ts_future;
  array[50] real ts_past;
  array[2] real y_init;            // initial measured population
  array[N, 2] real<lower = 0> y;   // measured population at measurement times
}

parameters {
  array[4] real<lower=0> theta; // { alpha, beta, gamma, delta }
  array[2] real<lower=0> z_init; // initial population
  array[2] real<lower=0> sigma; // measurement errors
}

transformed parameters {
  
  array[N, 2] real z
   = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-6, 1e-5, 1e3);


  array[50, 2] real z_future
    = integrate_ode_rk45(dz_dt, z[N], ts[N], ts_future,
                         theta, rep_array(0.0, 0), rep_array(0, 0),
                         1e-6, 1e-5, 1e3);
  array[50, 2] real z_past
    = integrate_ode_rk45(dz_dt, z_init, -50, ts_past,
                         theta, rep_array(0.0, 0), rep_array(0, 0),
                         1e-6, 1e-5, 1e3);
}

model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);                      // prior for initial population
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);   // likelihood for initial population
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);        // likelihood for observed population
  }
}

generated quantities {
  array[50, 2] real y_future_rep;
  array[50, 2] real y_past_rep;
  
  for (k in 1:2) {
    for (n in 1:50) {
      y_future_rep[n, k] = lognormal_rng(log(z_future[n, k]), sigma[k]);
      y_past_rep[n, k] = lognormal_rng(log(z_past[n, k]), sigma[k]);
    }
  }
}
