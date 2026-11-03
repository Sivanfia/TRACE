#!/usr/bin/env Rscript
# Only missing packages are installed. Existing packages are not upgraded.
required<-c('MASS','Matrix','glmnet','quadprog','Rcpp','RcppEigen','RcppArmadillo','spikeslab')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
print(data.frame(package=required,version=vapply(required,function(x)as.character(packageVersion(x)),'')))
cat('Full fitting requires a C++ compiler. Saved-result figures use base R only.\n')
