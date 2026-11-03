# Run once on the server in your chosen R environment. Requires CRAN access.
required <- c("MASS", "Matrix", "glmnet", "Rcpp", "RcppArmadillo", "RcppEigen", "quadprog", "spikeslab")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Installation incomplete: ", paste(missing, collapse = ", "))
print(data.frame(package = required, version = vapply(required, function(p) as.character(packageVersion(p)), "")))
cat("Rcpp compilation of TRADER and CONCERT requires the server's C++ toolchain; some R dependencies also require Fortran.\n")
cat("Headless Linux plotting requires R capabilities png=TRUE and cairo=TRUE.\n")
print(capabilities(c("png", "cairo")))
