# Exact original manuscript configuration, including rho_z=0.15, rho_c=0.3.
.real_root <- dirname(dirname(dirname(normalizePath(sys.frame(1)$ofile))))
CONFIG <- readRDS(file.path(.real_root,'retained_results/original_manifest.rds'))$config
rm(.real_root)
