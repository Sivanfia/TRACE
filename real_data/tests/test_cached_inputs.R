a<-commandArgs(TRUE);root<-normalizePath(a[1]);source(file.path(root,'vendor/trace_acvb_kernels.R'));source(file.path(root,'R/preprocess.R'))
source(file.path(root,'R/retained_inputs.R'))
files<-list.files(file.path(root,'retained_results/jobs'),pattern='\\.rds$',full.names=TRUE)
# A saved input must reconstruct the actual old split, not sample a new split.
for(f in files[seq(1,length(files),length.out=6)]){
 b<-readRDS(f);p<-reconstruct_real_inputs(b,file.path(root,'data/processed/gtex'))
 stopifnot(identical(p$split,b$split),identical(p$training$B,b$B),identical(p$precision,b$precision),is.null(p$training$beta_true),is.null(p$training$X_test),is.null(p$training$y_test))
 stopifnot(max(abs(colMeans(p$training$X)))<1e-10,max(abs(apply(p$training$X,2,sd)-1))<1e-10)
 # Target bootstrap, gene order, sample order and calibration are fixed.
 stopifnot(identical(colnames(p$training$X),b$feature_map$gene_id),identical(p$training$target_donor,b$split$train_donors))
 stopifnot(isTRUE(all.equal(p$test$y_log2_tpm,b$predictions$observed[b$predictions$method=='ACVB'],tolerance=1e-14)))
 bad<-b;bad$predictions$observed[bad$predictions$method=='ACVB']<-bad$predictions$observed[bad$predictions$method=='ACVB']+1
 stopifnot(inherits(try(reconstruct_real_inputs(bad,file.path(root,'data/processed/gtex')),silent=TRUE),'try-error'))
 cat('PASS: exact stored split and observed response alignment:',b$target,'rep',b$rep,'\n')
}
manifest<-readRDS(file.path(root,'retained_results/original_manifest.rds'))
check_retained_data_files(file.path(root,'data/processed/gtex'),manifest,manifest$config$targets)
tmp<-tempfile('changed_gtex_');dir.create(tmp);t<-manifest$config$targets[1];f<-paste0(t,'.cns_tpm.tsv.gz')
file.copy(file.path(root,'data/processed/gtex',f),file.path(tmp,f));con<-file(file.path(tmp,f),'ab');writeBin(as.raw(0),con);close(con)
stopifnot(inherits(try(check_retained_data_files(tmp,manifest,t),silent=TRUE),'try-error'))
unlink(tmp,recursive=TRUE)
cat('PASS: original data fingerprint enforced, including held-out predictor values\n')
