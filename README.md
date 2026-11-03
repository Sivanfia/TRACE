# TRACE

**Adaptive Credibility Reweighting for Robust Bayesian Model Selection with Source Summary Statistics**

TRACE combines individual-level target data with fitted coefficients from multiple source studies. It estimates source weights and uses separate latent indicators for target variable selection and coordinate-specific information transfer. Computation uses an approximated collapsed variational algorithm with variational EM updates for the source weights.

This repository contains the source code for the accompanying manuscript by Pengfei Chen, Yuehan Yang, Klinger Lu, and Chong You. It includes four simulation settings, a six-tissue GTEx application, plotting scripts, tests, and environment records. The statistical source files are identical to those in the complete reproducibility package, including the corrected zero-inclusion residual calculation.

## Source code and complete materials

The `main` branch contains source files. Large inputs, saved fits, results, and generated figures are distributed separately as the Release asset **`simulation_and_real_data_analysis.zip`**. Download it from this repository's **Releases** page. GitHub's automatically generated “Source code” ZIP contains only this source tree; it does not include that asset.

For the saved-result workflow below, extract the complete package into a **separate directory**, then restore the required inputs and result tables into your source checkout:

```bash
cd /path/to/TRACE
Rscript scripts/restore_release_data.R /path/to/simulation_and_real_data_analysis
```

Use the extracted directory containing `MANIFEST.json`, `simulation/`, and `real_data/` as the argument. The restoration script verifies the archived source version, copies the required external files, and leaves source files unchanged. It also stops before copying if an existing data file differs. Repeating restoration with identical files is safe.

The restored inputs and results are excluded from Git by `.gitignore`. Most archived simulation checkpoints and duplicate plots are unnecessary for the supplied workflows and are left in the complete Release. Alternatively, the complete extracted package can be run directly using its own README.

## Reproduce figures and tables

After restoration, run:

```bash
bash run_analysis.sh figures
```

This rebuilds displays from the saved replicate-level results without fitting models. It checks the numerical summaries against the archived manuscript data and writes new files to `outputs/reproduced/`.

Base R with PDF and PNG graphics support is sufficient. Key outputs include:

- `figures/setting01_main_f1.pdf`: Setting 1 at SNR 4, with six panels.
- `figures/setting02_main_f1.pdf`: both TRADER benchmark designs, including CONCERT.
- `figures/setting03_main_f1.pdf`: Setting 3 at SNR 3, with six panels.
- `figures/setting04_main_f1.pdf`: Setting 4 at noise standard deviation 3.
- `figures/real_data_all_tissues_mse.pdf`: the six-tissue comparison.
- `tables/`: MSE, model-size, and supplementary tables in CSV and TeX format.
- `checks/result_counts_and_agreement.csv`: result counts and agreement checks.

These paths are relative to `outputs/reproduced/`. Simulation error bars are one Monte Carlo standard error. The real-data MSE tables report standard deviations across splits.

## Refit the analyses

Full fitting requires R, a C++ compiler, and MASS, Matrix, glmnet, quadprog, Rcpp, RcppEigen, RcppArmadillo, and spikeslab. Install missing packages with:

```bash
Rscript scripts/install_dependencies.R
```

The original experiments used R 4.4.1 on Linux; package versions and session information are in `environment/`. The installer does not pin versions or upgrade installed packages. Results of fresh fits can vary across R, package, compiler, and BLAS versions.

After restoring the external data, refit all methods on a 48-core server:

```bash
nohup bash run_48cores.sh all > run_all.log 2>&1 &
```

The four simulation settings run sequentially, followed by the real-data analysis. At most 48 replication/split workers run at once, with one BLAS/OpenMP thread per worker. Full simulation reruns regenerate the data and spike-and-slab source summaries. New fits go into `runs/`, separately from the archived results.

To use fewer workers:

```bash
CORES=8 bash run_analysis.sh all
```

To refit TRACE-EM and Oracle TRACE in simulation, and TRACE-EM in real data, while retaining the saved competitors:

```bash
CORES=48 bash run_analysis.sh trace
```

This route checks regenerated inputs and initialization against archived references before merging results. Use the original experiment environment for it. Real data have no Oracle TRACE comparison because true source weights are unknown.

To run only Setting 3 with all methods:

```bash
SIM_FRESH_INPUTS=TRUE SIM_NCORES=48 SIM_OUTPUT_BASE="$PWD/runs/simulation_all" \
  bash simulation/run_setting.sh 03 all
```

The publication plotting workflow requires 50 replications. Smaller fitting runs are available for inspection. Use a new output directory when changing code, inputs, or parameters, so incompatible checkpoints are not mixed.

## Experiments

| Setting | Design | Parameter cells | Replications per cell |
|---|---|---:|---:|
| 1 | Shared contamination of active source coefficients | 48 | 50 |
| 2 | TRADER Designs I and III | 20 | 50 |
| 3 | Coordinate-dependent locally correct sources | 48 | 50 |
| 4 | Complementary source combinations and shared contamination | 21 | 50 |

Settings 1, 3, and 4 use 500 observations per source. Setting 2 uses 200 in Design I and 500 in Design III. There are 6,850 simulation replications and 41,100 method-level records, covering TRACE-EM, Oracle TRACE, ACVB, Trans-Lasso, TRADER, and CONCERT. Setting 4 uses contamination fractions from 0 to 0.6, correlation 0.3, and noise standard deviation 3.

The GTEx analysis uses six brain tissues, 50 donor-level splits per tissue, and five methods. Target test donors are excluded from every source. Target training donors can occur in source tissues. Selected models are refitted by target-training ridge regression with inner cross-validation for prediction evaluation. Data acquisition and preprocessing provenance are under `real_data/provenance/`.

## Directory map

```text
simulation/
  01_TRACE_identifiable_mixture_contamination/   # Setting 1
  02_TRACE_support_swap_contamination/           # Setting 2
  03_TRACE_mosaic_contamination/                 # Setting 3
  04_TRACE_complementary_sources_fixed_noise/    # Setting 4
  _common/                                     # kernels and method adapters
  tests/
real_data/
  R/, vendor/, config/                         # TRACE rerun with retained competitors
  full_analysis/                               # complete five-method analysis
  provenance/, tests/
manuscript/figures/                             # plotting source only
scripts/                                       # restoration, reconstruction, validation
environment/                                   # original environment records
provenance/                                    # source manifest and archive notes
run_analysis.sh                                # common entry point
run_48cores.sh                                 # 48-worker entry point
```

Data directories appear after restoration. Manuscript snapshots restored from the Release document the numerical reference used by these scripts; they are not a replacement for the current submission manuscript.

## Validation

The restoration helper can be tested without downloading the Release:

```bash
Rscript scripts/test_restore_release_data.R
```

After restoration, run the original validation suite:

```bash
bash run_analysis.sh validate
```

It checks source syntax, record counts, paired comparisons, saved real-data checkpoints, TRACE identities, input reconstruction, and result merging. A broader check performs small temporary fits:

```bash
FULL_VALIDATION=1 bash run_analysis.sh validate
```

See `VALIDATION.md` for checks performed when preparing this source distribution. Historical checks for the complete archive are preserved separately in `provenance/full_package_validation.md`. Neither packaging check reruns all formal 50-replication experiments.

## Attribution and licensing

Original comments, third-party notices, upstream identifiers, and adaptation notes are preserved; see `NOTICE.md`, `simulation/CONCERT_BASELINE.md`, and `simulation/TRANS_LASSO_OFFICIAL.md`. No open-source license is specified here for the authors' original code. Third-party code and data remain subject to their original terms.

If using these materials in research, please cite the accompanying TRACE manuscript and the relevant baseline methods. The manuscript is not assigned a publication DOI in this repository.
