# TRACE

**Adaptive Credibility Reweighting for Robust Bayesian Model Selection with Source Summary Statistics**

TRACE combines individual-level target data with fitted coefficients from multiple source studies. It learns source weights and uses separate latent indicators for target variable selection and coordinate-specific information transfer. Computation uses an approximated collapsed variational algorithm with variational EM updates for source weights.

This project provides the code and reproducibility materials for the accompanying manuscript by Pengfei Chen, Yuehan Yang, Klinger Lu, and Chong You. The complete package contains four simulation settings and a six-tissue GTEx gene-expression application.

## Download the complete package

Download `simulation_and_real_data_analysis.zip` from the **Releases** section of this repository. The archive includes source code, analysis inputs, saved results, plotting scripts, and environment records.

The commands below should be run from the root of the extracted package:

```bash
unzip simulation_and_real_data_analysis.zip
cd simulation_and_real_data_analysis
```

## Reproduce the manuscript figures and tables

To rebuild the figures and tables from the saved replication-level results:

```bash
bash run_analysis.sh figures
```

This mode requires **R with PDF and PNG graphics support**. It recomputes summaries, checks agreement with the manuscript data, and writes the outputs to `outputs/reproduced/`. It does not refit the statistical methods or overwrite the archived results.

The output folders contain:

- `figures/`: the four main simulation figures, supplementary simulation figures, transfer diagnostics, and GTEx comparisons.
- `tables/`: GTEx prediction results and supplementary simulation tables.
- `data/`: regenerated numerical summaries.
- `checks/`: result counts and agreement checks.

Simulation error bars represent one Monte Carlo standard error. The GTEx MSE table reports the mean and standard deviation across donor splits.

## Refit the analyses

The reported fits used **R 4.4.1 on Linux**. Original package versions and session information are recorded in `environment/`.

Full fitting requires Bash, R, a C++ compiler, and the following R packages:

`MASS`, `Matrix`, `glmnet`, `quadprog`, `Rcpp`, `RcppEigen`, `RcppArmadillo`, and `spikeslab`.

Install missing packages with:

```bash
Rscript scripts/install_dependencies.R
```

The installer adds missing packages without upgrading existing ones; it does not recreate the recorded package versions automatically.

To run all simulation methods and the full real-data comparison on a machine with 48 available CPU cores:

```bash
nohup bash run_48cores.sh all > run_all.log 2>&1 &
```

The four simulation settings run sequentially, followed by the GTEx application. The runner uses up to 48 workers and limits each worker to one BLAS/OpenMP thread. New fits are stored under `runs/`, separately from the archived results.

For fewer workers, use:

```bash
CORES=8 bash run_analysis.sh all
```

To refit only TRACE while retaining the saved competitor results:

```bash
CORES=48 bash run_analysis.sh trace
```

This refits TRACE-EM and Oracle TRACE in simulation, and TRACE-EM in the GTEx application. Oracle TRACE is used only in simulation, where the true source weights are known. This mode checks reconstructed inputs before combining new fits with the retained results and should use the recorded original environment.

Fresh fits may differ across R, package, compiler, and numerical-library versions. The saved-result reproduction mode provides the numerical results reported in the manuscript.

## Experiments

| Setting | Design | Parameter combinations | Replications per combination |
|---|---|---:|---:|
| 1 | Shared contamination on active coordinates | 48 | 50 |
| 2 | TRADER benchmark Designs I and III | 20 | 50 |
| 3 | Coordinate-dependent source compatibility | 48 | 50 |
| 4 | Complementary sources with shared contamination | 21 | 50 |

Simulation comparisons include TRACE-EM, Oracle TRACE, ACVB, Trans-Lasso, TRADER, and CONCERT. Settings 1, 3, and 4 use 500 observations per source. Setting 2 uses 200 per source in Design I and 500 in Design III. The package contains 6,850 simulation replications and 41,100 method-level records, including CONCERT in both TRADER designs.

The GTEx application uses six brain tissues with 50 donor-level train/test splits per target and five methods. Target test donors are excluded from all sources. Selected models are refitted by target-training ridge regression with an internally cross-validated penalty for prediction evaluation.

## Package structure

```text
simulation/        Simulation designs, method implementations, and saved results
real_data/         GTEx inputs, preprocessing, method implementations, and results
scripts/           Dependency installation, figure reconstruction, and validation
environment/       Original R sessions and package versions
provenance/        Data sources, implementation provenance, and result mapping
manuscript/        Manuscript snapshots, numerical summaries, and plotting sources
outputs/           Reproduced figures, tables, and checks
run_analysis.sh    Main entry point: figures, validate, trace, or all
run_48cores.sh     Runner for full analysis with 48 workers
```

Keep this directory structure when running the complete package. Plot reconstruction uses the files in `manuscript/data/` and `manuscript/figures/` as well as the saved analysis results.

## Validation

Run the standard checks with:

```bash
bash run_analysis.sh validate
```

For additional checks that perform small temporary fits:

```bash
FULL_VALIDATION=1 bash run_analysis.sh validate
```

The scope and results of the package assembly checks are documented in `VALIDATION.md`. They include syntax checks, result pairing and counts, numerical agreement, TRACE matrix identities, and input reconstruction. These checks do not rerun the complete 50-replication analyses.

## Attribution and licensing

Third-party implementations and GTEx data provenance are documented in `NOTICE.md` and the associated `provenance/` folders. Original attribution and notices are retained.

A license for the original TRACE code has not yet been specified. Third-party code and data remain subject to their original terms.

## Citation

Please cite the accompanying manuscript when reporting work that uses TRACE:

Chen, P., Yang, Y., Lu, K., and You, C. *TRACE: Adaptive Credibility Reweighting for Robust Bayesian Model Selection with Source Summary Statistics*. Manuscript.

