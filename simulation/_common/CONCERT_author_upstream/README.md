# CONCERT: Bayesian Transfer Learning

This repository contains the R code for the `Concert` algorithm, simulations, real-data applications for the linear and logistic models.
The paper, Covariate-Elaborated Robust Partial Information Transfer with Conditional Spike-and-Slab Prior, can be found at https://arxiv.org/abs/2404.03764.

- **Linear**
  - The main functions are in `Concert.R` .
  - Simulation examples are provided in `simu`.
  - Application on the **GTEx** data is provided in `realdata`.
- **Binary (logistic)**
  - The main functions are in `Concert-binary.R` .
  - Simulation examples are provided in `simu`.
  - Application on the **Lending Club** data is provided in `realdata`.

- The key function `Concert()` provides the following results:
  - `m_beta0`: Estimated slab mean for $\beta_0$.
  - `gamma_Z`: Estimated variational posterior inclusion probability for the target.
  - `m_betak`: Estimated slab mean for $\beta_k$.
  - `gamma_I`: Estimated variational posterior transferable probability for each source.
  - `V_beta0`: Estimated slab variance for $\beta_0$.
  - `V_betak`: Estimated slab variance for $\beta_k$.

