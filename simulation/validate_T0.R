args <- commandArgs(TRUE)
entry <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root <- if (length(args)) args[1] else dirname(normalizePath(entry))

read_defs <- function(path) {
  ans <- list()
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]])) {
      ans[[as.character(expr[[2]])]] <- expr
    }
  }
  ans
}

load_kernels <- function(path, core_path=NULL) {
  e <- new.env(parent=globalenv())
  load_defs <- function(defs) {
    for (expr in defs) {
      rhs <- expr[[3]]
      if (is.call(rhs) && identical(rhs[[1]],as.name('function'))) eval(expr,e)
      else if (is.symbol(rhs) && exists(as.character(rhs),e,inherits=TRUE)) eval(expr,e)
    }
  }
  if (!is.null(core_path)) load_defs(read_defs(core_path))
  load_defs(read_defs(path))
  e
}

# This oracle evaluates the complete matrices and the reduced model directly.
# Reintroducing the removed-coordinate source-prior term must fail this test.
direct_residual <- function(w,g,mu,dat,tau0,remove=integer(0),free_j=NULL) {
  keep <- setdiff(seq_along(w),remove)
  if (!length(keep)) return(list(R=dat$yTy,ld=0))
  tx <- get_ns_eff_for_oracle(dat)*dat$x_norm2/dat$n
  precision <- g*tx+(1-g)*tau0
  Om <- w %o% w
  diag(Om) <- w
  M <- dat$XtX*Om+diag(precision+1e-10,length(w))
  if (!is.null(free_j)) M[free_j,free_j] <- M[free_j,free_j]-1e-10
  rhs <- w*dat$Xty+g*tx*mu
  Mk <- M[keep,keep,drop=FALSE]
  list(R=as.numeric(dat$yTy+sum(g[keep]*tx[keep]*mu[keep]^2)-
       crossprod(rhs[keep],solve(Mk,rhs[keep]))),
       ld=as.numeric(determinant(Mk,logarithm=TRUE)$modulus))
}
get_ns_eff_for_oracle <- function(dat) dat$trace_source_precision_ns_eff

oracle_omega <- function(j,w,g,mu,dat,tau0,Ay,By,prior) {
  w1 <- w; w1[j] <- 1
  a <- direct_residual(w,g,mu,dat,tau0,remove=j)
  b <- direct_residual(w1,g,mu,dat,tau0,free_j=j)
  tx <- get_ns_eff_for_oracle(dat)*dat$x_norm2/dat$n
  prec <- g[j]*tx[j]+(1-g[j])*tau0
  lin <- prior-.5*(b$ld-a$ld-log(prec))+
    (Ay+dat$n/2)*log((By+a$R/2)/(By+b$R/2))
  pmin(pmax(plogis(lin),.Machine$double.eps),1-.Machine$double.eps)
}
oracle_gamma <- function(j,w,g,mu,dat,tau0,Ay,By,prior) {
  g0 <- g1 <- g; g0[j] <- 0; g1[j] <- 1
  a <- direct_residual(w,g0,mu,dat,tau0,free_j=j)
  b <- direct_residual(w,g1,mu,dat,tau0,free_j=j)
  tx <- get_ns_eff_for_oracle(dat)*dat$x_norm2/dat$n
  lin <- prior+.5*log(tx[j]/tau0)-.5*(b$ld-a$ld)+
    (Ay+dat$n/2)*log((By+a$R/2)/(By+b$R/2))
  pmin(pmax(plogis(lin),.Machine$double.eps),1-.Machine$double.eps)
}

set.seed(271027)
n <- 10L; p <- 4L
X <- matrix(rnorm(n*p),n,p); y <- rnorm(n)
dat <- list(X=X,y=y,n=n,p=p,K=2L,XtX=crossprod(X),Xty=as.numeric(crossprod(X,y)),
            x_norm2=colSums(X^2),yTy=sum(y^2),
            B=matrix(c(2,-1,1,.5,-1,2,-.5,1),p,2),
            trace_source_precision_ns_eff=40)
w <- c(.35,.6,.8,.45); mu <- as.numeric(dat$B%*%c(.65,.35))
tau0 <- 10; Ay <- By <- 1; pz <- -.8; pc <- -.3
files <- list.files(root,pattern='functions[.]R$',full.names=TRUE,recursive=TRUE)
files <- files[vapply(files,function(f) '.update_omega_all_woodbury' %in% names(read_defs(f)),logical(1))]
real_kernel <- file.path(root,'TRACE_real_data_with_CONCERT_server/vendor/trace_acvb_kernels.R')
if(file.exists(real_kernel))files <- c(files,real_kernel)
stopifnot(length(files)>0L)
failures <- character(0); checks <- 0L
assert_near <- function(got,want,label,tol=2e-7) {
  checks <<- checks+1L
  err <- max(abs(got-want))
  if (!is.finite(err)||err>tol) failures <<- c(failures,sprintf('%s: max error %.6g',label,err))
}
for (path in files) {
  core <- file.path(dirname(dirname(path)),'_common/borrow_functions_core.R')
  e <- load_kernels(path,if(file.exists(core))core else NULL)
  for (g in list(c(.2,.45,.7,.85),rep(0,p),rep(1,p))) {
    want <- vapply(1:p,oracle_omega,numeric(1),w=w,g=g,mu=mu,dat=dat,tau0=tau0,Ay=Ay,By=By,prior=pz)
    inp <- list(omega=w,gamma=g,bbar=mu,dat=dat,S0bar_inv=tau0,A_y=Ay,B_y=By,E_logit_rho_z=pz)
    for (nm in c('.update_omega_all_woodbury','.update_omega_all_cached_inverse','.update_omega_all_cached')) {
      if (exists(nm,e,inherits=FALSE)) assert_near(do.call(e[[nm]],inp),want,paste(path,nm))
    }
    got <- vapply(1:p,function(j)do.call(e$update_omega_j_old,c(list(j=j),inp)),numeric(1))
    assert_near(got,want,paste(path,'scalar'))
    seq_w <- w
    for (j in 1:p) seq_w[j] <- oracle_omega(j,seq_w,g,mu,dat,tau0,Ay,By,pz)
    if (exists('.update_omega_all_gauss_seidel',e,inherits=FALSE)) {
      assert_near(do.call(e$.update_omega_all_gauss_seidel,inp),seq_w,paste(path,'sequential'))
    }
    if (exists('.update_gamma_all_gauss_seidel',e,inherits=FALSE)) {
      seq_g <- g
      for (j in 1:p) seq_g[j] <- oracle_gamma(j,w,seq_g,mu,dat,tau0,Ay,By,pc)
      ginp <- inp; ginp$E_logit_rho_z <- NULL; ginp$E_logit_rho_c <- pc
      assert_near(do.call(e$.update_gamma_all_gauss_seidel,ginp),seq_g,paste(path,'sequential gamma cache'))
    }
    seq_w <- w; seq_g <- g
    for (j in 1:p) {
      seq_w[j] <- oracle_omega(j,seq_w,seq_g,mu,dat,tau0,Ay,By,pz)
      seq_g[j] <- oracle_gamma(j,seq_w,seq_g,mu,dat,tau0,Ay,By,pc)
    }
    if (exists('.update_omega_gamma_interleaved_gauss_seidel',e,inherits=FALSE)) {
      ans <- do.call(e$.update_omega_gamma_interleaved_gauss_seidel,c(inp,list(E_logit_rho_c=pc)))
      assert_near(ans$omega,seq_w,paste(path,'interleaved omega'))
      assert_near(ans$gamma,seq_g,paste(path,'interleaved gamma'))
    }
    # Inspect T1 from the actual vectorized body and compare it to a direct solve.
    env <- list2env(inp,parent=e)
    invisible(eval(body(e$.update_omega_all_woodbury),env))
    want1 <- vapply(1:p,function(j){ww<-w;ww[j]<-1;direct_residual(ww,g,mu,dat,tau0,free_j=j)$R},numeric(1))
    assert_near(env$T1,want1,paste(path,'included residual T1'),tol=1e-7)
  }
}
cat('Files:',length(files),'Checks:',checks,'Failures:',length(failures),'\n')
if (length(failures)) {
  cat(paste(head(failures,12),collapse='\n'),'\n')
  stop('TRACE zero-inclusion residual differs from direct reduced-model oracle.',call.=FALSE)
}
cat('PASS: reduced T0, unchanged included T1, all available inclusion backends and borrowing endpoints.\n')
