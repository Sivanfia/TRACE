require(here)

load(here("Linear","realdata","GTEx","GTEx.RData"))

gene.tar <- c("SH2D2A")
gene.tar.id <- which(colnames(data_all)==gene.tar)
ids <- setdiff(1:ncol(data_all),c(gene.tar.id,ncol(data_all)))
corrs <- sapply(ids, function(k){
  cor(data_all[,k],data_all[,gene.tar.id])
})
unselect.id <- ids[which(abs(corrs)<0.25)]
data_final <- data_all[,setdiff(1:ncol(data_all),unselect.id)]

Tissue.tar <- c("brain_amygdala", "brain_anterior_cingulate_cortex_ba24",
                "brain_cortex", "brain_hippocampus",
                "brain_nucleus_accumbens_basal_ganglia", "brain_substantia_nigra")

###################### load methods ##########################
require(MASS)
require(caret)
require(dplyr)
require(glmnet)
require(foreach)
require(doParallel)

source(here("Linear", "Concert.R"))
source(here("Linear", "compare", "TransLasso-main", "TransLasso-functions.R"))
require(glmtrans)
require(sparsevb)

for (j in 1:6){
  tissue.tar <- Tissue.tar[j]
  print(tissue.tar)
  gene.tar.id <- which(colnames(data_final)==gene.tar)
  tissue.tar.id <- which(tissue.names==tissue.tar)
  target.sample.id <- which(data_final$tissue==tissue.tar)
  n.tar <- length(target.sample.id)
  
  q_k <- 0.1; q_0 <- 0.1
  seed_num <- 11
  if (j==3){seed_num <- 999}
  if (j==4){seed_num <- 555}
  if (j==6){seed_num <- 1234}
  
  res <- c()
  for (kk in 1:10){
    print(kk)
    numCores <- 5
    registerDoParallel(numCores)
    results.kk <- foreach(i=((kk-1)*5+1):(kk*5), .packages = c("sparsevb", "glmnet", "glmtrans")) %dopar%{
        if (j==1){set.seed(seed_num*(i+24))}
        if (j==2){set.seed(seed_num*(i+2))}
        if (j==3){set.seed(seed_num*(i+10))}
        if (j==4){set.seed(seed_num*(i+32))}
        if (j==5){set.seed(seed_num*(i+28))}
        if (j==6){set.seed(seed_num*(i+44))}
        
        train.id <- sample(target.sample.id,size = 0.8 * n.tar)
        test.id <- setdiff(target.sample.id, train.id)
        gene.cons.id.train <- which(apply(data_final[train.id,1:(ncol(data_final)-1)],2,sd)==0)
        if (length(gene.cons.id.train)==0){
          data_new <- data_final
        }else{
          data_new <- data_final[,-gene.cons.id.train]
        }
        data_train <- as.matrix(rbind(data_new[train.id,1:(ncol(data_new)-1)],
                                      data_new[-target.sample.id,1:(ncol(data_new)-1)]))
        data_test <- as.matrix(data_new[test.id,1:(ncol(data_new)-1)])
        
        X <- data_train[,-gene.tar.id]
        y <- data_train[,gene.tar.id]
        X_test <- data_test[,-gene.tar.id]
        y_test <- data_test[,gene.tar.id]
        
        n_vec<- c(length(train.id),tissue.n[-tissue.tar.id])
        p <- ncol(X)
        n0 <- n_vec[1]; K <- length(n_vec)-1
        
        #################### Spike-and-slab parameters ####################
        tau=rep(1,K); eta=1
        
        start <- Sys.time()
        lasso_0 <- cv.glmnet(X[1:n0,], y[1:n0], family="gaussian", alpha=0)
        lasso_beta_0 <- as.numeric(lasso_0$glmnet.fit$beta[,which(lasso_0$lambda==lasso_0$lambda.min)])
        q_k <- 1
        for (k in 2:(K+1)){
          lasso_k <- cv.glmnet(X[sample_ind(n_vec, k)[[k]],], y[sample_ind(n_vec, k)[[k]]], family="gaussian", alpha=0)
          lasso_beta_k <- as.numeric(lasso_k$glmnet.fit$beta[,which(lasso_k$lambda==lasso_k$lambda.min)])
          q_k <- min(q_k, sum(abs(lasso_beta_0-lasso_beta_k) <= 0.008)/p)
        }

        Concert_res <- Concert(X, y, n_vec, tau=tau, q_k=q_k, eta=eta, q_0=q_0)
        end <- Sys.time()
        diff_time <- difftime(end, start, units="secs")[[1]]
        m_beta0 <- Concert_res[[1]]
        gamma_Z <- Concert_res[[2]]
        m_beta0 <- m_beta0*gamma_Z
        PMSE_Concert <- unlist(eval_pre(m_beta0, X_test, y_test))
        cat("Trial ", i, ": ", PMSE_Concert, " with time ", diff_time, " secs.", "\n", sep="")
        
        ###########Sparse VB#############
        Target_res_l <- svb.fit(X[1:n0,], y[1:n0], family="linear", slab="laplace")
        beta0_target_l <- Target_res_l$mu
        beta0_target_l <- beta0_target_l*Target_res_l$gamma
        PMSE_target_l <- unlist(eval_pre(beta0_target_l, X_test, y_test))
        
        Target_res <- solo(X[1:n0,], y[1:n0], eta=eta, q_0=q_0)
        beta0_target <- Target_res$m_beta0
        beta0_target <- beta0_target*Target_res$gamma_Z
        PMSE_target <- unlist(eval_pre(beta0_target, X_test, y_test))
        
        ###########Lasso#############
        res.lasso <- cv.glmnet(X[1:n0,], y[1:n0], family="gaussian")
        beta0_lasso <- as.numeric(res.lasso$glmnet.fit$beta[,which(res.lasso$lambda==res.lasso$lambda.min)])
        PMSE_lasso <- unlist(eval_pre(beta0_lasso, X_test, y_test))
        
        ###########Trans-Lasso#############
        l1=T
        prop.re1 <- Trans.lasso(X, y, n_vec, I.til = 1:50, l1 = l1)
        prop.re2 <- Trans.lasso(X, y, n_vec, I.til = 101:n_vec[1], l1=l1)
        beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
        PMSE_TL <- unlist(eval_pre(beta.prop, X_test, y_test))
        
        ###########GLMtrans#############
        target <- list(x=X[1:n0,], y=y[1:n0])
        source <- list()
        for (k in 2:(K+1)){
          source_k <- list(x=X[sample_ind(n_vec, k)[[k]],], y=y[sample_ind(n_vec, k)[[k]]])
          source[[k-1]] <- source_k
        }
        glmtrans_fit <- glmtrans(target, source, intercept = F)
        beta_glmtrans <- glmtrans_fit$beta[-1]
        PMSE_glmtrans <- unlist(eval_pre(beta_glmtrans, X_test, y_test))
        
        list(PMSE_Concert=PMSE_Concert, PMSE_target_l=PMSE_target_l, PMSE_target = PMSE_target,
             PMSE_lasso=PMSE_lasso, PMSE_TL=PMSE_TL, PMSE_glmtrans=PMSE_glmtrans)
    }
    stopImplicitCluster()
    res <- c(res, results.kk)
    measure_mat <- matrix(unlist(res),ncol=6,byrow = TRUE)
    colnames(measure_mat) <- c("Concert", "SparseVB", "NaiveVB", "Lasso", "TransLasso", "TransGLM")
    ntrial <- nrow(measure_mat)
    measure_mat <- rbind(measure_mat,apply(measure_mat[1:ntrial,], 2, mean))
    measure_mat <- rbind(measure_mat,apply(measure_mat[1:ntrial,], 2, sd))
    print(measure_mat[(ntrial+1),])
    
    OUTPUT_DIR <- here("Linear","realdata","res")
    if (!dir.exists(OUTPUT_DIR)) {
      dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }
    loc <- paste(c(OUTPUT_DIR, "/", tissue.tar, ".rda"), collapse ="")
    save_name <- paste(c(tissue.tar, " seed ", seed_num), collapse = "")
    save_res <- list(save_name, measure_mat)
    save(save_res, file=loc)
  }
}


for (i in 1:6){
  print(i)
  tissue.tar <- Tissue.tar[i]
  OUTPUT_DIR <- here("Linear","realdata","res")
  loc <- paste(c(OUTPUT_DIR, "/", tissue.tar,".rda"), collapse ="")
  load(loc)
  res_mat <- save_res[[2]]
  ntrial <- nrow(res_mat)-2
  print(round(apply(res_mat[1:ntrial,], 2, mean),3))
  print(round(apply(res_mat[1:ntrial,], 2, sd),3))
}


