functions{
  array[] real dz_dt(real t, // time
                     array[] real z,
                     // system state {prey, predator}
                     array[] real theta, // parameters
                     array[] real x_r, // unused data
                     array[] int x_i) {
    real u = z[1];
    real v = z[2];
    
    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];
    
    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return {du_dt, dv_dt};
}
}

data {
  int<lower=0> N; // number of measurement times
  array[N] real ts; // measurement times > 0
}

transformed data{
    // ALEATORIOS independientes DE LA PREVIA (1), (2), (3)
    // (1) parametros del sistema LV
    real alpha_sim = normal_rng(1, 0.5);
    real gamma_sim = normal_rng(1, 0.5);
    real beta_sim = normal_rng(0.05, 0.05);
    real delta_sim = normal_rng(0.05, 0.05);
    // (2) escala verosimilitudes
    array[2] real<lower=0> sigma_sim = lognormal_rng(rep_vector(-1,2), 1);
    // (3) condiciones iniciales
    array[2] real<lower=0> z_init_sim = lognormal_rng(rep_vector(log(10),2), 1);
    
    // Almacena los valores generados de pars LV en vector theta_sim
    array[4] real<lower=0> theta_sim;
    theta_sim[1] = alpha_sim;
    theta_sim[2] = beta_sim;
    theta_sim[3] = gamma_sim;
    theta_sim[4] = delta_sim;

    // enchufa las simulaciones de theta arriba en la edo y resuelve
    // para 20 años. Genera sinteticos usando las soluciones
    // para simular de la posterior de theta (vector de parámetros).
    // Estos datos se pasarán para simular de la posterior
    array[N, 2] real z_sim = integrate_ode_rk45(dz_dt, z_init_sim, 0, ts, theta_sim,
                                          rep_array(0.0, 0), rep_array(0, 0), 
                                          1e-6, 1e-5, 1e7);

    // simula datos sinteticos (todo esto equivaldría a obtener otros
    // datos, e.g que la compañia recolectora de pieles hubiese hecho esto más veces.
    // siempre que se siga cumpliendo este modelo.)
    array[2] real<lower=0>  y_init;
    array[N, 2] real<lower=0>  y;
    for (k in 1 : 2) {
    y_init[k] = lognormal_rng(log(z_init_sim[k]), sigma_sim[k]);
    y[ : , k] = lognormal_rng(log(z_sim[ : , k]), sigma_sim[k]);
    }
}

parameters {
  array[4] real<lower=0> theta; // { alpha, beta, gamma, delta }
  array[2] real<lower=0> z_init; // initial population
  array[2] real<lower=0> sigma; // measurement errors
}

// ¿necesito lo siguiente?
// sí, son funciones phi_n = f(a,b,c,d,z_init) necesarias para el
// modelo de verosimilitud
transformed parameters {
  array[N, 2] real z = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                                          rep_array(0.0, 0), rep_array(0, 0), 
                                          1e-6, 1e-5, 1e7);
}
// aunque solo vayamos a usar muestras posteriores de theta, 
// el modelo se mantiene (pues la posterior se ve influenciada 
// también por las y's i.e la verosimilitud)
model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
  for (k in 1 : 2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
    y[ : , k] ~ lognormal(log(z[ : , k]), sigma[k]);
  }
}
generated quantities {
    int<lower=0, upper=1> alpha_lt_sim = theta[1] < theta_sim[1];
    int<lower=0, upper=1> beta_lt_sim = theta[2] < theta_sim[2];
    int<lower=0, upper=1> gamma_lt_sim = theta[3] < theta_sim[3];
    int<lower=0, upper=1> delta_lt_sim = theta[4] < theta_sim[4];
}
