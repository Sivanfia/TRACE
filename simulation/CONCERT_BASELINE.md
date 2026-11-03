# CONCERT baseline integration

## Backend and provenance

The bundle vendors the authors' linear CONCERT implementation from
`RuqianZhang/CONCERT`, pinned at commit
`fea55703c39655a196d3aa5296d1851f0fd53547`. The upstream `Concert.cpp` is
compiled with the algebraically equivalent optimizations documented in
`_common/CONCERT_author_upstream/OPTIMIZATION_NOTE.md`. The local adapter validates dimensions,
records the commit, and converts the returned posterior quantities into the
simulation package's common result schema.

The upstream source and third-party notice are under
`_common/CONCERT_author_upstream/`.

## Inputs

CONCERT is fit to the stacked target and raw source data:

```text
X_all = rbind(X_target, X_source[[1]], ..., X_source[[K]])
y_all = c(y_target, y_source[[1]], ..., y_source[[K]])
n_vec_trans = c(n_target, n_source_1, ..., n_source_K)
```

It does not use TRACE's coefficient-summary matrix `B` as a substitute for
individual-level source observations. This matches the CONCERT likelihood and
the way Trans-Lasso already uses raw source observations in this bundle.

Settings 01/03/04 default to 500 observations per source. TRACE and TRADER
receive the same original-scale spikeslab/gnet summaries fitted from these observations;
CONCERT and Trans-Lasso receive the corresponding raw observations. The
coefficient summaries preserve a different amount of information from the
raw data, but their underlying observations and sample counts are shared.

## Default hyperparameters

- `tau = 1` for every source.
- `eta = 1`.
- `q_0 = 0.1`.
- CAVI threshold `1e-6`, maximum 1000 iterations.
- Target support is selected by `gamma_Z >= 0.5`.
- The reported coefficient is the official posterior mean
  `m_beta0 * gamma_Z`.

The default `q_k` is the authors' simulation rule. Ridge CV is run separately
on the target and every source; source `k` agreement is the fraction of
coordinates satisfying `abs(beta_ridge_target - beta_ridge_source_k) <= 0.008`.
The common `q_k` is the minimum agreement over sources. A numerical clamp to
`[1e-6, 1-1e-6]` only prevents infinite logits when the empirical fraction is
exactly zero or one; both the raw and used values are saved.

For a sensitivity run, set `concert_qk_rule = "fixed"` and choose
`concert_qk_fixed` in each setting's `config.R`.

## Metrics and diagnostics

Prediction MSE uses the same independent observed-response test data as all
other methods. F1, precision, recall, TP, FP, FN, and model size use the common
0.5 inclusion threshold. Runtime includes the `q_k` ridge calibration and the
CONCERT fit.

The bundled optimized-equivalent backend returns the iteration count,
convergence flag, and final coefficient-update difference. The adapter saves
these actual diagnostics; it does not manufacture an ELBO value.

Each replicate saves `concert_diagnostics.csv` and
`concert_coordinate_posterior.csv`. The run manifest hashes the vendored C++
source as well as the R adapter.
