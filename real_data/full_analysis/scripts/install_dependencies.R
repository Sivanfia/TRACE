packages<-c('glmnet','Rcpp','RcppArmadillo','RcppEigen','spikeslab','quadprog')
missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
if(any(!vapply(packages,requireNamespace,logical(1),quietly=TRUE)))stop('Dependency installation incomplete')
print(vapply(packages,function(p)as.character(packageVersion(p)),''))
