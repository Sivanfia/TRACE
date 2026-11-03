# Exact configuration of results_50rep, which generated the manuscript table.
# rho_z=0.15, rho_c=0.3; target bootstrap and refit CV settings are unchanged.
# Read the saved R object to preserve every numeric value without text rounding.
.config_root <- dirname(dirname(normalizePath(sys.frame(1)$ofile)))
CONFIG <- readRDS(file.path(.config_root,'retained_results/original_manifest.rds'))$config
rm(.config_root)
