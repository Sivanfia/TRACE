# Official Trans-Lasso baseline

The four settings use the authors' Q-aggregation implementation at commit
`285009edf3e638567474108f823adcc1dc6c7c4a`.  The pinned source is
`TransLasso-author-285009e.R` in each setting and is loaded through
`translasso_official.R`.

The adapter uses deterministic two-fold target aggregation, preserves the
author algorithm, and stops on backend errors; it does not silently switch to
the former ranked-source fallback.  Performance tables record backend,
commit, source hash, adapter hash, split, numerical policy, and RNG policy.

Trans-Lasso and CONCERT use raw source-level observations.  TRACE and the
default TRADER run use the coefficient-summary matrix `B`; see
`CONCERT_BASELINE.md` for the information-budget qualification.
