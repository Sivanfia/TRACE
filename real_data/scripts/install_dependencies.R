packages<-c('glmnet','quadprog')
missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
cat('Dependencies ready. Existing packages were not upgraded.\n')
