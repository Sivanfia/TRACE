# Third-party methods and data

This package preserves the method source files and provenance supplied with the original analyses. Attribution, upstream commit identifiers and data acquisition notes remain with those files. No new license is asserted for third-party code or datasets.

- Trans-Lasso: author implementation and pinned source notes in `simulation/TRANS_LASSO_OFFICIAL.md`, the setting-level `TransLasso-author-285009e.R` files, and `real_data/full_analysis/provenance/translasso_commit.json`.
- CONCERT: pinned upstream source and the optimized equivalent linear backend in `simulation/_common/CONCERT_author_upstream/`; adaptation notes in `simulation/CONCERT_BASELINE.md`. The real-data adapter and upstream files are under `real_data/full_analysis/vendor/`.
- TRADER: the R and C++ bridge/core files supplied with each setting and under `real_data/full_analysis/vendor/`.
- ACVB and TRACE: the kernels and analysis wrappers supplied with this project, including the zero-inclusion-residual correction used for the reported results.
- GTEx v8 tissue-module input files: provenance under `real_data/provenance/data_provenance/` and `real_data/full_analysis/provenance/data_provenance/`.

Original source attribution and comments have been retained. Package-level edits concern portable entry points, complete result integration, output labels and validation; they do not change the statistical settings used for the saved results.
