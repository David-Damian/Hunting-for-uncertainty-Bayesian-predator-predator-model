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
generated quantities {
    array[4] real<lower=0> theta; // { alpha, beta, gamma, delta }
    array[2] real z_init; // initial population
    array[2] real sigma; // measurement errors
    array[2] real y_init_rep;
    array[2] real y_rep;
    array[1] real ts;

    theta[1] = normal_rng(1, 0.5);
    theta[2] = normal_rng(0.05, 0.05);
    theta[3] = normal_rng(1, 0.5);
    theta[4] = normal_rng(0.05, 0.05);

    z_init[1] = lognormal_rng(log(10), 1);
    z_init[2] = lognormal_rng(log(10), 1);

    sigma[1] = lognormal_rng(-1, 1);
    sigma[2] = lognormal_rng(-1, 1);

    ts[1] = 1;
    

    array[1, 2] real z = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                                            rep_array(0.0, 0), 
                                            rep_array(0, 0), 1e-5, 1e-3, 5e2);
    for (k in 1 : 2) {
        y_init_rep[k] = lognormal_rng(log(z_init[k]), sigma[k]);
        y_rep[k] = lognormal_rng(log(z[1, k]), sigma[k]);
    }
}
