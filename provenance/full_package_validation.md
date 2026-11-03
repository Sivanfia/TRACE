# Reproducibility validation

Assembly date: 28 October 2026. Validation environment: R 4.6.0 on macOS. The reported computations used R 4.4.1 on Linux; original versions and session information are preserved in `environment/`.

Passed checks:

- All 113 R source files in the execution/validation directories parse; root shell entry points pass shell syntax checks.
- Four complete simulation grids: 48, 20, 48 and 21 cells, with 50 paired replications and six methods per cell (41,100 method-level records).
- The recovered 500 CONCERT Design III records have unique keys; 1,000 matched TRADER/Trans-Lasso records confirm the corresponding historical seeds and performance metrics.
- Simulation means, standard deviations, Monte Carlo standard errors and effective counts recomputed from replication tables agree with the current manuscript summaries (numerical tolerance 1e-9).
- All 1,500 real-data metric records agree with the 300 saved five-method checkpoints. The reproduced MSE and model-size tables agree with the manuscript data.
- Supplementary Setting 2 metric ranges and Setting 4 diagnostics regenerated from the saved results agree with the active appendix table entries.
- Current manuscript snapshots are byte-identical to the supplied local main text and supplement; their comments are retained.
- Corrected TRACE kernels are byte-identical to the kernels used for the corrected results, including the full real-data runner's kernel.
- Matrix identity checks pass for the reduced zero-inclusion residual, included-state compensation and implemented inclusion/transfer backends.
- Stored GTEx splits, observations, feature order and preprocessing reconstruct correctly; input fingerprint and result-merging rejection checks pass.
- Temporary correction-only fitting tests pass for Settings 1, 3 and 4 and both TRADER designs, including resume and retained-comparator checks.
- All six simulation methods fit the same temporary small dataset. The full real-data entry reconstructs raw inputs and fits TRACE. Six temporary tissue fits preserve the saved competitor rows and produce tables and figures.
- The reported plots and tables regenerate successfully from the package. Main Setting 2 and 4 figures and the six-tissue real-data figure were visually inspected.
- Saved-result reproduction also passes after moving the entry points to a directory with spaces in its path.

The full 50-replication fitting analyses were not rerun during assembly. The saved reported results are included. Temporary small-fitting tests validate execution and numerical identities; they are not substituted for publication results. To rerun the complete analyses, use `bash run_analysis.sh all` under the recorded environment. To rebuild the reported figures and tables, use `bash run_analysis.sh figures`.

Logs are in `provenance/validation_full.log`, `reproduce_results.log`, `relocation_check.log` and `result_integration.log`. Historical manifests and input fingerprints in the archived results describe the original fits; Design III recovery and presentation updates are documented separately in `provenance/assembly_changes.json`.
