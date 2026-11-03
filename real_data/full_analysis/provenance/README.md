# Code provenance

TRACE and ACVB numerical kernels were extracted from the user's delivered
`TRACE_four_settings_01_to_04_server/04_TRACE_complementary_sources_fixed_noise`.
`simulation_source_sha256.json` records the upstream file hashes. The extraction
removes evaluation fields requiring simulated truth; the application never
passes a fake beta_true, alpha_true, tag_true, or held-out outcome into a fit.
Initializers exposed by the application only allow non-oracle initialization.
The EM alpha and coordinate-update equations are unchanged.

TRADER uses that same package's residual O(n p) C++ sampler. The bridge's
truth-dependent selection metric call and output field were removed; MCMC,
source screening, selection rule and posterior coefficient mean are unchanged.

Trans-Lasso author repository: https://github.com/saili0103/TransLasso
File revision: 285009edf3e638567474108f823adcc1dc6c7c4a
Downloaded 2026-09-15. `translasso_commit.json` records the GitHub response;
`TransLasso-original.R` retains the downloaded source. Local modifications in
`TransLasso-functions.R`:

1. All-zero aggregation returns the same list shape as the nonzero case.
2. Subtract the maximum log weight before exponentiating (identical normalized
   weights mathematically, avoids floating-point underflow).
3. Bound marginal screening count by p (avoids NA indices for small p).
4. Set intercept=FALSE and standardize=FALSE for the input scale used by all
   five application methods; the application standardizes outer training data.
5. Omit the unused beta-MSE helper. No algorithmic fallback occurs implicitly.

CONCERT author repository: https://github.com/RuqianZhang/CONCERT
File revision: fea55703c39655a196d3aa5296d1851f0fd53547. The complete upstream
linear directory is retained in `vendor/CONCERT_author_upstream/Linear`, and
`concert_commit.json` records the commit and SHA-256 checksums. The upstream
`Concert.cpp` is compiled without algorithmic edits. `vendor/concert_bridge.R`
only supplies a portable path, reproduces the authors' GTEx ridge-agreement
rule for `q_k`, clamps exact boundary probabilities by 1e-6 to avoid undefined
log odds, validates the result, and maps `m_beta0 * gamma_Z` to the shared
prediction interface. The upstream commit did not contain a license file; see
the bundled third-party notice.

Data: GTEx v8 public tissue TPM, extracted using the CNS module supplied in
https://github.com/RuqianZhang/CONCERT . The original project's download manifest
and gene-mapping audit are retained in `data_provenance/`. The public author
GTEx.RData was not available; these inputs are rebuilt and are not an exact
reproduction of their historical gene-mapping version or preprocessing.

Source files are included for the user's research reproduction; preserve
upstream attribution and review upstream licensing before further redistribution.

Real-data compatibility fixes: ns retrieval uses exact `[["n_vec", exact=TRUE]]`
to prevent R partial matching to `n_vec_trans`; this preserves the intended
source precision. TRADER switches spikeslab's bigp.smalln option off when n>=p
(the package rejects that path for n>p), without changing the estimator.
ACVB detects a Jacobi two-cycle and then applies 0.5 relaxation, retaining the
same coordinate equations and post-selection refit. Convergence also requires
the undamped map residual < tolerance. The original no-relaxation implementation
is retained and is reproduced when acvb_adaptive_damping=FALSE; this equivalence
is checked in tests/validate_pipeline.R.
