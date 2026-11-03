# Official Trans-Lasso baseline

This setting uses the author implementation from
<https://github.com/saili0103/TransLasso>, pinned at commit
`285009edf3e638567474108f823adcc1dc6c7c4a`.

`TransLasso-author-285009e.R` is the preserved author source. The local
`translasso_official.R` adapter calls `Trans.lasso(..., l1=TRUE)` with the
author's candidate ranking, cross-validation and Q-aggregation. It performs
the planned two target folds and averages the two estimates. It fails loudly
if the pinned file is absent, modified or returns invalid coefficients; there
is no fallback estimator.

The adapter contains two numerical compatibility fixes: stable normalization
of Q-aggregation log weights and a valid list return when every candidate
coefficient is zero. These do not change the estimator in regular cases.
The adapter also restores the caller's RNG state after the author's internal
cross-validation.

For the current TRADER-paper DGP, Trans-Lasso receives the same raw target and
ten raw source samples generated for the other methods. The run manifest
stores its backend identity and source hashes. Output includes response
prediction MSE and support metrics; coefficient MSE is intentionally absent.
