# Code for demo results
## numerical results
## plots

########################################################
########################################################
################# numerical results ####################
########################################################
########################################################
require(glmnet)
require(foreach)
require(doParallel)
require(here)


source(here("Linear", "Concert.R"))

################## setting 2 ######################
Coef_gen_high <- function(p, K, ratio, sig_beta, q, sig_delta, exact=T){
  beta0 <- c(sig_beta, rep(0, p-length(sig_beta)))
  beta <- matrix(0, nrow=p, ncol=K)
  delta <- matrix(0, nrow=p, ncol=K)
  
  s <- length(sig_beta)
  for (k in 1:K){
    beta[sample(1:s, ceiling(s*ratio)),k] <- beta0[sample(1:s, ceiling(s*ratio))]
    I_delta_k <- sample((1+s):p, q, replace=F)
    if (exact){
      delta[I_delta_k, k] <- sig_delta*(1-2*rbinom(q,1,0.5))
    } else{
      delta[I_delta_k, k] <- rnorm(q, beta0[I_delta_k], 1.2)
    }
    beta[,k] <- beta[,k] + delta[,k]
  }
  return(list(beta0=beta0, beta=beta, delta=delta))
}

numCores <- max(detectCores()-3, 1)
registerDoParallel(numCores)

s_set <- c(4)
ratio_set <- c(0.5)

for (s in s_set){
  for (ratio in ratio_set){
    cat("s: ", s, " ratio: ", ratio, '\n', sep = '')
    set.seed(11+s*10+ratio*10)
    p <- 10
    n0 <- 150
    K <- 5
    sig_beta <- rep(0.5, s)
    q <- 3; sig_delta <- 1
    
    n_vec <- c(n0, rep(200, K))
    Sig_X <- diag(1, p)
    coefs <- Coef_gen_high(p, K, ratio, sig_beta, q, sig_delta, exact=T)
    beta00 <- coefs$beta0
    delta0 <- coefs$delta
    betak0 <- coefs$beta
    true_coefs <- coefs
    B0 <- cbind(beta00, coefs$beta)

    ntrial <- 100
    measure_mat <- matrix(0, nrow=ntrial+2, ncol=5)
    colnames(measure_mat) <- c("rmse_Concert","prmse_Concert","TPR_Concert","FDR_Concert","time_Concert")
    res_all <- list()
    loop_res <- foreach(trial=1:ntrial, .combine = rbind) %dopar% {
      #################### generate the data ####################
      set.seed(11+63*7+trial*7)
      
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
      
      #################### Spike-and-slab parameters ####################
      tau=rep(1,K); eta=1; q_0=0.1
      q_k <- 0.5
      
      start <- Sys.time()
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
      
      list(c(rmse_Concert,prmse_Concert,TPR_Concert,FDR_Concert,time_Concert),
           concert_res)
    }
    
    for (trial in 1:ntrial){
      measure_mat[trial,] <- loop_res[trial,][[1]]
      res_all[[trial]] <- loop_res[trial,][[2]]
    }
    measure_mat[ntrial+1,] <- apply(measure_mat[1:ntrial,], 2, mean)
    measure_mat[ntrial+2,] <- apply(measure_mat[1:ntrial,], 2, var)
    res_all[[ntrial+1]] <- B0
    
    OUTPUT_DIR <- here("Linear", "simu", "res", "demo")
    if (!dir.exists(OUTPUT_DIR)) {
      dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }
    write.csv(measure_mat, paste(OUTPUT_DIR, "/Res setting 2 demo s ", s, " ratio ", ratio," sig beta ", sig_beta[1], ".csv", sep=""), row.names = F) 
    save(res_all, file=paste(OUTPUT_DIR, "/Res setting 2 demo s ", s, " ratio ", ratio," sig beta ", sig_beta[1], ".rda", sep=""))
  }
}
stopImplicitCluster()


######################## setting 3 ########################
Coef_gen_high <- function(p, K, sig_beta, size_A0, h, sig_delta1, q, sig_delta2, exact=T){
  beta0 <- c(sig_beta, rep(0, p-length(sig_beta)))
  beta <- matrix(rep(beta0, K), nrow=p, ncol=K)
  delta <- matrix(0, nrow=p, ncol=K)
  for (k in 1:K){
    if (k <= size_A0){
      if (exact){
        # I_delta_k <- sample(1:p, h, replace=F) # difference loc
        I_delta_k <- sample((1+length(sig_beta)):p, h, replace=F)
        delta[I_delta_k, k] <- -sig_delta1*(1-2*rbinom(h,1,0.5))
      } else{
        delta[1:100, k] <- rnorm(100, 0, h/100)
      }
      beta[,k] <- beta[,k] + delta[,k]
    } else{
      if (exact){
        # I_delta_k <- sample(1:p, q, replace=F)
        I_delta_k <- sample((1+length(sig_beta)):p, q, replace=F)
        delta[I_delta_k, k] <- -sig_delta2*(1-2*rbinom(q,1,0.5))
      } else{
        delta[1:100, k] <- rnorm(100, 0, q/100)
      }
      beta[,k] <- beta[,k] + delta[,k]
    }
  }
  return(list(beta0=beta0, beta=beta, delta=delta))
}

numCores <- max(detectCores()-3, 1)
registerDoParallel(numCores)


q_set <- c(3)
sig_delta2_set <- c(1)
for (q in q_set){
  for (sig_delta2 in sig_delta2_set){
    cat("q: ", q, " sig_delta2: ", sig_delta2, '\n', sep = '')
    set.seed(11+q+sig_delta2*10)
    p <- 10 
    n0 <- 150
    K <- 5
    s <- 4
    size_A0 <- 0
    sig_beta <- rep(0.5, s)
    h <- 2; sig_delta1 <- 0.5
    
    n_vec <- c(n0, rep(200, K))
    Sig_X <- diag(1, p)
    coefs <- Coef_gen_high(p, K, sig_beta, size_A0, h, sig_delta1, q, sig_delta2,
                           exact=T)
    beta00 <- coefs$beta0
    delta0 <- coefs$delta
    betak0 <- coefs$beta
    true_coefs <- coefs
    B0 <- cbind(beta00, coefs$beta)
    
    ntrial <- 100
    measure_mat <- matrix(0, nrow=ntrial+2, ncol=5)
    colnames(measure_mat) <- c("rmse_Concert","prmse_Concert","TPR_Concert","FDR_Concert","time_Concert")
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
      
      #################### Spike-and-slab parameters ####################
      tau=rep(1,K); eta=1; q_0=0.1
      
      start <- Sys.time()
      q_k <- 0.5
      
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
      
      list(c(rmse_Concert,prmse_Concert,TPR_Concert,FDR_Concert,time_Concert),
           concert_res)
    }
    
    for (trial in 1:ntrial){
      measure_mat[trial,] <- loop_res[trial,][[1]]
      res_all[[trial]] <- loop_res[trial,][[2]]
    }
    measure_mat[ntrial+1,] <- apply(measure_mat[1:ntrial,], 2, mean)
    measure_mat[ntrial+2,] <- apply(measure_mat[1:ntrial,], 2, var)
    res_all[[ntrial+1]] <- B0
    
    OUTPUT_DIR <- here("Linear", "simu", "res", "demo")
    if (!dir.exists(OUTPUT_DIR)) {
      dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }
    write.csv(measure_mat, paste(OUTPUT_DIR, "/Res setting 3 demo q", q," sig_delta2 ", sig_delta2, ".csv", sep=""), row.names = F)  
    save(res_all, file=paste(OUTPUT_DIR, "/Res setting 3 demo q", q," sig_delta2 ", sig_delta2, ".rda", sep=""))
  }
}
stopImplicitCluster()




########################################################
########################################################
####################### plots ##########################
########################################################
########################################################
library(ggplot2)
library(patchwork)


#################### Demo setting 2 ####################
p <- 10
K <- 5
s <- c(4)
ratio <- c(0.5)
params_all <- as.matrix(expand.grid(s, ratio))
result_setting2 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  s <- params_all[ii,1]
  ratio <- params_all[ii,2]
  params <- c(s, ratio)
  
  loc <- here("Linear", "simu", "res", "demo")
  load(paste(loc, "/Res setting 2 demo s ", s, " ratio ", ratio, " sig beta 0.5.rda", sep=""))
  true_para <- res_all[[101]]
  beta0_0 <- true_para[,1]
  betak_0 <- true_para[,2:6]
  Z <- (beta0_0!=0)*1
  I_mat <- (betak_0==0)*1
  for (k in 1:K){I_mat[1:s,k] <- (betak_0[1:s,k]==beta0_0[1:s])*1}
  
  probs <- matrix(0, nrow=(K+1)*p, ncol=102)
  for (trial in 1:100){
    trial_res <- res_all[[trial]]
    gamma_Z <- trial_res$gamma_Z
    gamma_I <- as.numeric(trial_res$gamma_I)
    probs[,trial] <- c(gamma_Z, gamma_I)
  }
  probs[,101] <- apply(probs[,1:100],1,mean)
  probs[,102] <- apply(probs[,1:100],1,sd)
  prob_res <- cbind(probs[,101], c(Z*2, as.numeric(I_mat)), rep(0:K, each=p), rep(1:p,(K+1)))
  prob_res <- data.frame(prob_res)
  colnames(prob_res) <- c("pip","positive","source","index")
  prob_res$positive <- factor(prob_res$positive,
                              levels=c(2,1,0), labels=c("T-Pos","S-Pos","Neg"))
  result_setting2 <- rbind(result_setting2, prob_res)
}

######## Plot ########
data_measure_setting2 <- as.data.frame(result_setting2)
source_labels <- c(paste("Source ", K:1, sep=""), "Target  ")
data_measure_setting2$source <- factor(data_measure_setting2$source,
                                           levels = K:0, labels = source_labels)
covariate_index <- 1:p
data_measure_setting2$index <- factor(data_measure_setting2$index,
                                       levels = 1:p, labels = paste("X",1:p,sep=""))

labels <- paste("X[", 1:p, "]", sep = "")


p1_1<-ggplot(data_measure_setting2, aes(x = index, y = source, fill = pip))+
  labs(x ="(b)", y = "",fill = "Prob")+
  geom_tile(width = 1, height = 0.95, color = "white", size = 0.2) + 
  scale_fill_gradientn(colours = c("skyblue2", "yellow", "brown2"),
                       values = scales::rescale(c(-0.5, -0.25, 0, 0.25, 0.5)))+
  scale_x_discrete(labels = parse(text = labels)) +
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.title = element_text(size=15),
        text = element_text(size = 22),
        axis.text.x=element_text(size = 20),
        panel.background = element_rect(fill = "white", colour = "white"),  # Ensure background is white
        plot.background = element_rect(fill = "white", colour = "white"),   # Ensure plot background is white
        legend.key.height=unit(0.6, 'cm'),
        legend.key.width=unit(1.2, 'cm'))

p2_1<-ggplot(data_measure_setting2, aes(x = index, y = source, fill = positive))+
  labs(x ="(a)", y =  expression(bold('Case 1:')~' Heterogeneous'),fill = "Identity")+
  geom_tile(width = 1, height = 0.95, color = "white", size = 0.2) +
  scale_fill_manual(values=c("brown2","#ff7010","skyblue2"))+
  scale_x_discrete(labels = parse(text = labels)) +
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 20),
        axis.text.y=element_text(size = 20),
        panel.background = element_rect(fill = "white", colour = "white"),  # Ensure background is white
        plot.background = element_rect(fill = "white", colour = "white"),   # Ensure plot background is white
        legend.title = element_text(size=20),
        text = element_text(size = 22))

p3 <- p2_1+p1_1


#################### Demo setting 3 ####################
p <- 10
K <- 5
q <- c(3)
sig_delta2 <- c(1)
s <- c(4)
params_all <- as.matrix(expand.grid(q, sig_delta2))
result_setting3 <- c()
#### combine results ####
for (ii in 1:nrow(params_all)){
  q <- params_all[ii,1]
  sig_delta2 <- params_all[ii,2]
  params <- c(q, sig_delta2)
  loc <- here("Linear", "simu", "res", "demo")
  load(paste(loc, "/Res setting 3 demo q", q, " sig_delta2 ", sig_delta2, ".rda", sep=""))
  true_para <- res_all[[101]]
  beta0_0 <- true_para[,1]
  betak_0 <- true_para[,2:6]
  Z <- (beta0_0!=0)*1
  I_mat <- (betak_0==0)*1
  for (k in 1:K){I_mat[1:s,k] <- (betak_0[1:s,k]==beta0_0[1:s])*1}
  
  probs <- matrix(0, nrow=(K+1)*p, ncol=102)
  for (trial in 1:100){
    trial_res <- res_all[[trial]]
    gamma_Z <- trial_res$gamma_Z
    gamma_I <- as.numeric(trial_res$gamma_I)
    probs[,trial] <- c(gamma_Z, gamma_I)
  }
  probs[,101] <- apply(probs[,1:100],1,mean)
  probs[,102] <- apply(probs[,1:100],1,sd)
  prob_res <- cbind(probs[,101], c(Z*2, as.numeric(I_mat)), rep(0:K, each=p), rep(1:p,(K+1)))
  prob_res <- data.frame(prob_res)
  colnames(prob_res) <- c("pip","positive","source","index")
  prob_res$positive <- factor(prob_res$positive,
                              levels=c(2,1,0), labels=c("Target Positive","Source Positive","Negative"))
  
  result_setting3 <- rbind(result_setting3, prob_res)
}

######## Plot ########
data_measure_setting3 <- as.data.frame(result_setting3)

source_labels <- c(paste("Source ", K:1, sep=""), "Target  ")
data_measure_setting3$source <- factor(data_measure_setting3$source,
                                       levels = K:0, labels = source_labels)
covariate_index <- 1:p
data_measure_setting3$index <- factor(data_measure_setting3$index,
                                      levels = 1:p, labels = paste("X",1:p,sep=""))

labels <- paste("X[", 1:p, "]", sep = "")

p1<-ggplot(data_measure_setting3, aes(x = index, y = source, fill = pip))+
  labs(x ="(d)", y ="",fill = "Prob")+
  geom_tile(width = 1, height = 0.95, color = "white", size = 0.2) + 
  scale_fill_gradientn(colours = c("skyblue2", "yellow", "brown2"),
                       values = scales::rescale(c(-0.5, -0.25, 0, 0.25, 0.5)))+
  scale_x_discrete(labels = parse(text = labels)) +
  theme(legend.position = "bottom",
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size = 20),
        legend.title = element_text(size=20),
        text = element_text(size = 22),
        panel.background = element_rect(fill = "white", colour = "white"),  # Ensure background is white
        plot.background = element_rect(fill = "white", colour = "white"),   # Ensure plot background is white
        legend.key.height=unit(0.6, 'cm'),
        legend.key.width=unit(1.2, 'cm'))

p2<-ggplot(data_measure_setting3, aes(x = index, y = source, fill = positive))+
  labs(x ="(c)", y =  expression(bold('Case 2:')~' Redundant'),fill = "Identity")+
  geom_tile(width = 1, height = 0.95, color = "white", size = 0.2) + 
  scale_fill_manual(values=c("brown2","#ff7010","skyblue2"))+
  scale_x_discrete(labels = parse(text = labels)) +
  theme(legend.position = "bottom",
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 20),
        axis.text.y=element_text(size = 20),
        panel.background = element_rect(fill = "white", colour = "white"),  # Ensure background is white
        plot.background = element_rect(fill = "white", colour = "white"),   # Ensure plot background is white
        legend.title = element_text(size=20),
        text = element_text(size = 22))

p4 <- p2+p1

p5 <- p3/p4


OUTPUT_DIR <- here("Linear", "simu", "pics")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
ggsave(p5, filename=paste0(OUTPUT_DIR,"/Demo-both.png",collapse=""), width=16, height=11)


