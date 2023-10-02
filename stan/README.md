# Stan models

This directory contains the [Stan](https://mc-stan.org/) models which we will use for the simulation experiments. All the models use the same likelihood function

$$
L(\theta|t, \eta, \beta) = \prod^{I}_{i = 1}[f(t_i)]^{\delta_O}[1 - F(t_i)]^{\delta_R}[F(t_i + t_{start}) - F(t_i)]^{\delta_I},
$$

where $\delta_O$, $\delta_R$, and $\delta_I$ are indicator variables for whether an observation is fully observed, right censored, or interval censored respectively. The function $f(.)$ is the Weibull PDF,

$$
f(t) = \beta \lambda t^{1 - \beta} \exp{[- \lambda t^\beta]},
$$

and the function $F(.)$ is the Weibull CDF,

$$
F(t) = 1 - \exp{[-\lambda t ^ \beta]}.
$$

The difference between the three models are the priors:

- Non-informative independent priors: `non-informative-model.stan`

- Independent informative priors: `independent-prior-model.stan`

- Joint informative Prior: `joint-prior-model.stan`

---

Stan Development Team. 2023. Stan Modelling Language Users Guide and Reference Manual, 2.21.0. <https://mc-stan.org>
