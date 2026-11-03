# Code for setting 1 for linear regression
## numerical results
## plots

########################################################
########################################################
################# numerical results ####################
########################################################
########################################################
require(glmnet)
require(ncvreg)
require(foreach)
require(doParallel)
require(here)

source(here("Linear", "Concert.R"))
source(here("Linear", "compare", "TransLasso-main", "TransLasso-functions.R"))
require(glmtrans)
require(sparsevb)

Coef_gen_high <- function(p, K, sig_beta, size_A0, h, sig_delta1, q, sig_delta2){
  beta0 <- c(sig_beta, rep(0, p-length(sig_beta))) # add intercept
  beta <- matrix(rep(beta0, K), nrow=p, ncol=K)
  delta <- matrix(0, nrow=p, ncol=K)
  for (k in 1:K){
    if (k <= size_A0){
      I_delta_k <- sample(1:p, h, replace=F) # difference loc
      delta[I_delta_k, k] <- -sig_delta1*(1-2*rbinom(h,1,0.5))
      beta[,k] <- beta[,k] + delta[,k]
    } else{
      I_delta_k <- sample(1:p, q, replace=F)
      delta[I_delta_k, k] <- -sig_delta2*(1-2*rbinom(q,1,0.5))
      beta[,k] <- beta[,k] + delta[,k]
    }
  }
  return(list(beta0=beta0, beta=beta, delta=delta))
}

numCores <- max(detectCores()-3, 1)
registerDoParallel(numCores)

h_set <- c(4,8,12)
size_A0_set <- c(0,2,4,6,8,10)
for (h in h_set){
for (size_A0 in size_A0_set){
  cat("h: ", h, " size_A0: ", size_A0, '\n', sep = '')
  set.seed(11+h*10+size_A0)
  p <- 200
  n0 <- 150
  K <- 10
  s <- 16
  sig_beta <- rep(0.5, s)
  sig_delta1 <- 0.5
  q <- 20; sig_delta2 <- 1
  
  n_vec <- c(n0, rep(100, K))
  Sig_X <- diag(1, p)
  coefs <- Coef_gen_high(p, K, sig_beta, size_A0, h, sig_delta1, q, sig_delta2)
  beta00 <- coefs$beta0
  delta0 <- coefs$delta
  betak0 <- coefs$beta
  true_coefs <- coefs
  B0 <- cbind(beta00, coefs$beta)
  
  ntrial <- 100
  measure_mat <- matrix(0, nrow=ntrial+2, ncol=30)
  colnames(measure_mat) <- c("rmse_Concert","rmse_Target_l","rmse_solo","rmse_Lasso","rmse_TL","rmse_GlmTrans",
                             "prmse_Concert","prmse_Target_l","prmse_solo","prmse_Lasso","prmse_TL","prmse_GlmTrans",
                             "TPR_Concert","TPR_Target_l","TPR_solo","TPR_Lasso","TPR_TL","TPR_GlmTrans",
                             "FDR_Concert","FDR_Target_l","FDR_solo","FDR_Lasso","FDR_TL","FDR_GlmTrans",
                             "time_Concert","time_Target_l","time_solo","time_Lasso","time_TL","time_GlmTrans")
  res_all <- list()
  loop_res <- foreach(trial=1:ntrial, .combine = rbind) %dopar% {
    #################### generate the data ####################
    set.seed(11+trial)
    X <- NULL
    y <- NULL
    for (k in 1:(K+1)){
      X <- rbind(X, MASS::mvrnorm(n_vec[k], rep(0, p), Sig_X))
      ind <- sample_ind(n_vec, k)[[k]]
      y <- c(y, X[ind,]%*%B0[,k] + rnorm(n_vec[k],0,1))
    }
    n_test <- n_vec[1]
    X_test <- MASS::mvrnorm(n_test, rep(0, p), Sig_X)
    y_test <- X_test%*%B0[,1] + rnorm(n_test,0,1)

    #################### Concert ####################
    tau=rep(1,K); eta=1; q_0=0.1
    
    start <- Sys.time()
    lasso_0 <- cv.glmnet(X[1:n0,], y[1:n0], family="gaussian", alpha=0)
    lasso_beta_0 <- as.numeric(lasso_0$glmnet.fit$beta[,which(lasso_0$lambda==lasso_0$lambda.min)])
    q_k <- 1
    for (k in 2:(K+1)){
      lasso_k <- cv.glmnet(X[sample_ind(n_vec, k)[[k]],], y[sample_ind(n_vec, k)[[k]]], family="gaussian", alpha=0)
      lasso_beta_k <- as.numeric(lasso_k$glmnet.fit$beta[,which(lasso_k$lambda==lasso_k$lambda.min)])
      q_k <- min(q_k, sum(abs(lasso_beta_0-lasso_beta_k) <= 0.008)/p)
    }

    concert_res <- Concert(X, y, n_vec, tau=tau, q_k=q_k, eta=eta, q_0=q_0)
    end <- Sys.time()
    time_Concert <- difftime(end, start, units="secs")[[1]]
    m_beta0 <- concert_res[[1]]
    gamma_Z <- concert_res[[2]]
    m_beta0 <- m_beta0*gamma_Z
    Z <- (gamma_Z>=0.5)*1
    rmse_Concert <- sqrt(sum((beta00-m_beta0)^2))
    TPR_Concert <- perf(Z, beta00)[[1]]
    FDR_Concert <- perf(Z, beta00)[[2]]
    prmse_Concert <- unlist(eval_pre(m_beta0, X_test, y_test))
    
    ########### Sparse VB #############
    start <- Sys.time()
    Target_res_l <- svb.fit(X[1:n0,], y[1:n0], family="linear", slab="laplace")
    end <- Sys.time()
    time_target_l <- difftime(end, start, units="secs")[[1]]
    beta0_target_l <- Target_res_l$mu
    Z_target_l <- (Target_res_l$gamma>=0.5)*1
    beta0_target_l <- beta0_target_l*Target_res_l$gamma
    rmse_target_l <- sqrt(sum((beta00-beta0_target_l)^2))
    TPR_target_l <- perf(Z_target_l, beta00)[[1]]
    FDR_target_l <- perf(Z_target_l, beta00)[[2]]
    prmse_target_l <- unlist(eval_pre(beta0_target_l, X_test, y_test))
    
    ########### solo #############
    eta=1; q_0=0.1
    start <- Sys.time()
    Target_res <- solo(X[1:n0,], y[1:n0], eta=eta, q_0=q_0)
    end <- Sys.time()
    time_target <- difftime(end, start, units="secs")[[1]]
    beta0_target <- Target_res$m_beta0
    Z_target <- (Target_res$gamma_Z>=0.5)*1
    beta0_target <- beta0_target*Target_res$gamma_Z
    rmse_target <- sqrt(sum((beta00-beta0_target)^2))
    TPR_target <- perf(Z_target, beta00)[[1]]
    FDR_target <- perf(Z_target, beta00)[[2]]
    prmse_target <- unlist(eval_pre(beta0_target, X_test, y_test))

    ########### Lasso #############
    start <- Sys.time()
    res.lasso <- cv.glmnet(X[1:n0,], y[1:n0], family="gaussian")
    end <- Sys.time()
    time_lasso <- difftime(end, start, units="secs")[[1]]
    beta0_lasso <- as.numeric(res.lasso$glmnet.fit$beta[,which(res.lasso$lambda==res.lasso$lambda.min)])
    rmse_lasso <- sqrt(sum((beta00-beta0_lasso)^2))
    TPR_lasso <- perf(beta0_lasso, beta00)[[1]]
    FDR_lasso <- perf(beta0_lasso, beta00)[[2]]
    prmse_lasso <- unlist(eval_pre(beta0_lasso, X_test, y_test))

    ########### Trans-Lasso #############
    start <- Sys.time()
    l1=T
    prop.re1 <- Trans.lasso(X, y, n_vec, I.til = 1:50, l1 = l1)
    prop.re2 <- Trans.lasso(X, y, n_vec, I.til = 101:n_vec[1], l1=l1)
    beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
    end <- Sys.time()
    time_TL <- difftime(end, start, units="secs")[[1]]
    rmse_TL <- sqrt(sum((beta00-beta.prop)^2))
    TPR_TL <- perf(beta.prop, beta00)[[1]]
    FDR_TL <- perf(beta.prop, beta00)[[2]]
    prmse_TL <- unlist(eval_pre(beta.prop, X_test, y_test))
    
    ########### TransGLM #############
    target <- list(x=X[1:n0,], y=y[1:n0])
    source <- list()
    for (k in 2:(K+1)){
      source_k <- list(x=X[sample_ind(n_vec, k)[[k]],], y=y[sample_ind(n_vec, k)[[k]]])
      source[[k-1]] <- source_k
    }
    start <- Sys.time()
    glmtrans_fit <- glmtrans(target, source, intercept = F, detection.info = F)
    end <- Sys.time()
    time_glmtrans <- difftime(end, start, units="secs")[[1]]
    beta_glmtrans <- glmtrans_fit$beta[-1]
    rmse_glmtrans <- sqrt(sum((beta00-beta_glmtrans)^2))
    TPR_glmtrans <- perf(beta_glmtrans, beta00)[[1]]
    FDR_glmtrans <- perf(beta_glmtrans, beta00)[[2]]
    prmse_glmtrans <- unlist(eval_pre(beta_glmtrans, X_test, y_test))

    list(c(rmse_Concert, rmse_target_l, rmse_target, rmse_lasso, rmse_TL, rmse_glmtrans, 
           prmse_Concert, prmse_target_l, prmse_target, prmse_lasso, prmse_TL, prmse_glmtrans,
           TPR_Concert, TPR_target_l, TPR_target, TPR_lasso, TPR_TL, TPR_glmtrans,
           FDR_Concert, FDR_target_l, FDR_target, FDR_lasso, FDR_TL, FDR_glmtrans,
           time_Concert, time_target_l, time_target, time_lasso, time_TL, time_glmtrans),
         concert_res)
  }
  
  for (trial in 1:ntrial){
    measure_mat[trial,] <- loop_res[trial,][[1]]
    res_all[[trial]] <- loop_res[trial,][[2]]
  }
  measure_mat[ntrial+1,] <- apply(measure_mat[1:ntrial,], 2, mean)
  measure_mat[ntrial+2,] <- apply(measure_mat[1:ntrial,], 2, var)
  res_all[[ntrial+1]] <- B0

  OUTPUT_DIR <- here("Linear", "simu", "res", "setting 1")
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(measure_mat, paste(OUTPUT_DIR, "/Res setting 1 h ", h, " size_A0 ", size_A0, ".csv", sep=""), row.names = F)
  save(res_all, file=paste(OUTPUT_DIR, "/Res setting 1 h ", h, " size_A0 ", size_A0, ".rda", sep=""))
}
}
stopImplicitCluster()


########################################################
########################################################
####################### plots ##########################
########################################################
########################################################
library(ggplot2)
library(tidyr)
library(Rmisc)
library(patchwork)


#################### estimation error ####################
h <- c(4,8,12)
A0 <- c(0,2,4,6,8,10)
Method <- c("Concert","SparseVB","NaiveVB", "Lasso", "TransLasso", "TransGLM")
params_all <- as.matrix(expand.grid(h, A0))
result_setting1 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  h <- params_all[ii,1]
  A <- params_all[ii,2]
  params <- c(h,A)
  
  loc <- here("Linear", "simu", "res", "setting 1")
  loc_res <- read.csv(paste(loc, "/Res setting 1 h ", h, " size_A0 ", A, ".csv", sep=""))
  loc_rmse <- as.numeric(loc_res[101,c(1:6)])
  loc_se <- as.numeric(loc_res[102,c(1:6)])
  loc_measure <- as.data.frame(cbind(A, loc_rmse, loc_se))
  h <- paste("h=",h,sep="")
  loc_measure <- cbind(h, loc_measure, Method)
  
  result_setting1 <- rbind(result_setting1, loc_measure)
}

data_measure_setting1 <- as.data.frame(result_setting1)

data_measure_setting1$h <- factor(data_measure_setting1$h,
                                  levels = c("h=4","h=8","h=12"),
                                  labels = c("Linear: h=4","Linear: h=8","Linear: h=12"))
data_measure_setting1$Method <- factor(data_measure_setting1$Method, levels = Method)

p1_Gaussian <-ggplot(data_measure_setting1, aes(x = A, y=loc_rmse, color = Method))+
  geom_line(aes(linetype=Method, color=Method), size = 0.5)+
  geom_point(aes(color=Method, shape=Method), size = 1)+
  geom_errorbar(aes(ymin=loc_rmse-loc_se,
                    ymax=loc_rmse+loc_se),
                width=0.5)+
  labs(x ="Informative set number |A|", y ="Estimation error")+
  scale_x_continuous(breaks = seq(0,10,2))+
  theme(plot.title = element_text(hjust = 0.5, size=10),
        text = element_text(size = 10),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),legend.position = "bottom",
        legend.key.size = unit(0.1, "inches"), legend.title = element_text(size=8))+
  facet_grid(~h)


OUTPUT_DIR <- here("Linear", "simu", "pics")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
ggsave(p1, filename=paste0(OUTPUT_DIR, "/setting", 1, "-est.png",collapse=""), width=4.2, height=2.4)


######################## prediction error #####################
h <- c(4,8,12)
A0 <- c(0,2,4,6,8,10)
Method <- c("Concert","SparseVB","NaiveVB", "Lasso", "TransLasso", "TransGLM")
params_all <- as.matrix(expand.grid(h, A0))
result_setting1 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  h <- params_all[ii,1]
  A <- params_all[ii,2]
  params <- c(h,A)
  
  loc <- here("Linear", "simu", "res", "setting 1")
  loc_res <- read.csv(paste(loc, "/Res setting 1 h ", h, " size_A0 ", A, ".csv", sep=""))
  loc_rmse <- as.numeric(loc_res[101,c(1:6)+6])
  loc_se <- as.numeric(loc_res[102,c(1:6)+6])
  loc_measure <- as.data.frame(cbind(A, loc_rmse, loc_se))
  h <- paste("h=",h,sep="")
  loc_measure <- cbind(h, loc_measure, Method)
  
  result_setting1 <- rbind(result_setting1, loc_measure)
}

data_measure_setting1 <- as.data.frame(result_setting1)

data_measure_setting1$h <- factor(data_measure_setting1$h,
                                  levels = c("h=4","h=8","h=12"),
                                  labels = c("Linear: h=4","Linear: h=8","Linear: h=12"))
data_measure_setting1$Method <- factor(data_measure_setting1$Method, levels = Method)

p1 <-ggplot(data_measure_setting1, aes(x = A, y=loc_rmse, color = Method))+
  geom_line(aes(linetype=Method, color=Method), size = 0.5)+
  geom_point(aes(color=Method, shape=Method), size = 1)+
  geom_errorbar(aes(ymin=loc_rmse-loc_se,
                    ymax=loc_rmse+loc_se),
                width=0.5)+
  labs(x ="Informative set number |A|", y ="Prediction error")+
  scale_x_continuous(breaks = seq(0,10,2))+
  theme(plot.title = element_text(hjust = 0.5, size=10),
        text = element_text(size = 10),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),legend.position = "bottom",
        legend.key.size = unit(0.1, "inches"), legend.title = element_text(size=8))+
  facet_grid(~h)


ggsave(p1, filename=paste0(OUTPUT_DIR, "/setting", 1, "-pred.png",collapse=""), width=4.2, height=2.4)



####################### tpr & fdr ##########################
##### Gaussian
h <- c(4,8,12)
A0 <- c(0,2,4,6,8,10)
Method <- c("Concert","SparseVB","NaiveVB", "Lasso", "TransLasso", "TransGLM")
params_all <- as.matrix(expand.grid(h, A0))
result_setting1 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  h <- params_all[ii,1]
  A <- params_all[ii,2]
  params <- c(h,A)
  
  loc <- here("Linear", "simu", "res", "setting 1")
  loc_res <- read.csv(paste(loc, "/Res setting 1 h ", h, " size_A0 ", A, ".csv", sep=""))
  loc_rmse <- as.numeric(loc_res[101,c(1:6)+12])
  loc_se <- as.numeric(loc_res[102,c(1:6)+12])
  loc_measure <- as.data.frame(cbind(A, loc_rmse, loc_se))
  h <- paste("h=",h,sep="")
  loc_measure <- cbind(h, loc_measure, Method)
  
  result_setting1 <- rbind(result_setting1, loc_measure)
}

data_measure_setting1 <- as.data.frame(result_setting1)

data_measure_setting1$h <- factor(data_measure_setting1$h,
                                  levels = c("h=4","h=8","h=12"),
                                  labels = c("Linear: h=4","Linear: h=8","Linear: h=12"))
data_measure_setting1$Method <- factor(data_measure_setting1$Method, levels = Method)

p1_Gaussian <-ggplot(data_measure_setting1, aes(x = A, y=loc_rmse, color = Method))+
  geom_line(aes(linetype=Method, color=Method), size = 0.5)+
  geom_point(aes(color=Method, shape=Method), size = 1)+
  geom_errorbar(aes(ymin=loc_rmse-loc_se,
                    ymax=loc_rmse+loc_se),
                width=0.5)+
  labs(x ="Informative set number |A|", y ="True positive rate")+
  scale_x_continuous(breaks = seq(0,10,2))+
  theme(plot.title = element_text(hjust = 0.5, size=10),
        text = element_text(size = 10),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),legend.position = "right",
        legend.key.size = unit(0.1, "inches"), legend.title = element_text(size=8))+
  facet_grid(~h)


h <- c(4,8,12)
A0 <- c(0,2,4,6,8,10)
Method <- c("Concert","SparseVB","NaiveVB", "Lasso", "TransLasso", "TransGLM")
params_all <- as.matrix(expand.grid(h, A0))
result_setting1 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  h <- params_all[ii,1]
  A <- params_all[ii,2]
  params <- c(h,A)
  
  loc <- here("Linear", "simu", "res", "setting 1")
  loc_res <- read.csv(paste(loc, "/Res setting 1 h ", h, " size_A0 ", A, ".csv", sep=""))
  loc_rmse <- as.numeric(loc_res[101,c(1:6)+18])
  loc_se <- as.numeric(loc_res[102,c(1:6)+18])
  loc_measure <- as.data.frame(cbind(A, loc_rmse, loc_se))
  h <- paste("h=",h,sep="")
  loc_measure <- cbind(h, loc_measure, Method)
  
  result_setting1 <- rbind(result_setting1, loc_measure)
}

data_measure_setting1 <- as.data.frame(result_setting1)

data_measure_setting1$h <- factor(data_measure_setting1$h,
                                  levels = c("h=4","h=8","h=12"),
                                  labels = c("Linear: h=4","Linear: h=8","Linear: h=12"))
data_measure_setting1$Method <- factor(data_measure_setting1$Method, levels = Method)

p1_Gaussian_fdr <-ggplot(data_measure_setting1, aes(x = A, y=loc_rmse, color = Method))+
  geom_line(aes(linetype=Method, color=Method), size = 0.5)+
  geom_point(aes(color=Method, shape=Method), size = 1)+
  geom_errorbar(aes(ymin=loc_rmse-loc_se,
                    ymax=loc_rmse+loc_se),
                width=0.5)+
  labs(x ="Informative set number |A|", y ="False discovery rate")+
  scale_x_continuous(breaks = seq(0,10,2))+
  theme(plot.title = element_text(hjust = 0.5, size=10),
        text = element_text(size = 10),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),legend.position = "right",
        legend.key.size = unit(0.1, "inches"), legend.title = element_text(size=8))+
  facet_grid(~h)


p1 <- p1_Gaussian/p1_Gaussian_fdr


ggsave(p1, filename=paste0(OUTPUT_DIR, "/setting", 1, "-tpr.pdf",collapse=""), width=5.6, height=3.6)




######################## similarity selection #############################
s <- 16
K <- 10
h <- c(4,8,12)
A0 <- c(0,2,4,6,8,10)
params_all <- as.matrix(expand.grid(h, A0))
result_setting1 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  h <- params_all[ii,1]
  A <- params_all[ii,2]
  params <- c(h,A)
  
  loc <- here("Linear", "simu", "res", "setting 1")
  load(paste(loc, "/Res setting 1 h ", h, " size_A0 ", A, ".rda", sep=""))
  true_para <- res_all[[101]]
  beta0_0 <- true_para[,1]
  betak_0 <- true_para[,2:11]
  I_mat <- (betak_0==0)*1
  for (k in 1:K){I_mat[1:s,k] <- (betak_0[1:s,k]==beta0_0[1:s])*1}
  
  probs <- matrix(0, nrow=200*10, ncol=102)
  for (trial in 1:100){
    trial_res <- res_all[[trial]]
    gamma_I <- as.numeric(trial_res$gamma_I)
    probs[,trial] <- gamma_I
  }
  probs[,101] <- apply(probs[,1:100],1,mean)
  probs[,102] <- apply(probs[,1:100],1,sd)
  prob_res <- cbind(1:2000,probs[,101], as.numeric(I_mat))
  prob_res <- data.frame(prob_res)
  colnames(prob_res) <- c("index","pip","positive")
  prob_res$positive <- factor(prob_res$positive)
  hlevel <- paste("h=",h,sep="")
  A0level <- paste("A0=",A,sep="")
  prob_res <- cbind(prob_res, hlevel, A0level)
  
  result_setting1 <- rbind(result_setting1, prob_res)
}

data_measure_setting1 <- as.data.frame(result_setting1)
data_measure_setting1$hlevel <- factor(data_measure_setting1$hlevel,
                                       levels = c("h=4","h=8","h=12"),
                                       labels = c("h=4","h=8","h=12"))
data_measure_setting1$A0level <- factor(data_measure_setting1$A0level,
                                        levels = c("A0=0","A0=2","A0=4","A0=6","A0=8","A0=10"),
                                        labels = c("|A|=0","|A|=2","|A|=4","|A|=6","|A|=8","|A|=10"))

p1 <- ggplot(data=data_measure_setting1,aes(x=index,y=pip,group=positive))+
  geom_point(aes(shape=positive,color=positive), size=0.3)+
  scale_color_manual(values=c("dodgerblue","orangered"))+
  labs(x="Covariate index",y="Posterior similarity probability",color="",shape="")+
  scale_shape_manual(values=c(3,1))+
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.position = "right",
        legend.key.size = unit(0.1, "inches"), 
        legend.title = element_text(size=10))+
  facet_grid(A0level~hlevel)


ggsave(p1, filename=paste0(OUTPUT_DIR, "/setting", 1, "-source.png",collapse=""), width=8, height=6)





