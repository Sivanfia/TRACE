## Target-only fixed-point solver: gamma is exactly zero and alpha is frozen.
## Start with the raw block map; after a detected two-cycle use eta=0.5.
## Convergence uses the undamped map residual at the returned omega.
trace_target_only_initialize <- function(dat, omega, S0bar_inv,
                                         rho_z, A_y=1, B_y=1,
                                         max_iter=1000L, tol=1e-6,
                                         verbose=FALSE) {
  stopifnot(length(omega)==dat$p, length(max_iter)==1L,is.finite(max_iter),
            max_iter>=1,length(tol)==1L,is.finite(tol),tol>0,
            is.finite(rho_z),rho_z>0,rho_z<1)
  max_iter <- as.integer(max_iter)
  rows <- vector("list",max_iter)
  converged <- FALSE
  initial_omega <- omega
  lag2 <- NULL
  damping_active <- FALSE
  damping_trigger_iter <- NA_integer_
  for(iter in seq_len(max_iter)) {
    proposed <- .update_omega_all_woodbury(
      omega=omega,gamma=rep(0,dat$p),bbar=rep(0,dat$p),dat=dat,
      S0bar_inv=S0bar_inv,A_y=A_y,B_y=B_y,E_logit_rho_z=qlogis(rho_z))
    if(!all(is.finite(proposed))) stop("Nonfinite target-only update.")
    residual <- max(abs(proposed-omega))
    cycle_distance <- if(is.null(lag2)) NA_real_ else max(abs(proposed-lag2))
    if(!damping_active && iter>=3L && residual>max(100*tol,1e-4) &&
       is.finite(cycle_distance) && cycle_distance<=max(1e-6,1e-3*residual)) {
      damping_active <- TRUE
      damping_trigger_iter <- iter
      if(verbose)cat(sprintf("[target only] detected two-cycle at iter=%d; activating eta=0.5\n",iter))
    }
    eta <- if(damping_active) .5 else 1
    # Keep the current point when it passes F(omega)-omega < tol.
    # This makes the saved residual a genuine residual of the returned omega.
    converged <- residual<tol
    omega_next <- if(converged) omega else omega+eta*(proposed-omega)
    step_delta <- max(abs(omega_next-omega))
    rows[[iter]] <- data.frame(iter=iter,max_omega_change=step_delta,
      raw_fixed_point_residual=residual,cycle_distance=cycle_distance,
      damping_active=damping_active,damping_factor=eta,
      size05=sum(omega_next>.5),sum_omega=sum(omega_next),gamma_max=0,alpha_updated=FALSE)
    lag2 <- omega
    omega <- omega_next
    if(verbose && (iter==1L || iter%%20L==0L || converged))
      cat(sprintf("[target only] iter=%d | raw residual=%.4e | size05=%d | eta=%.1f | gamma=0\n",iter,residual,sum(omega>.5),eta))
    if(converged)break
  }
  list(omega=omega,initial_omega=initial_omega,converged=converged,n_iter=iter,
       max_change=residual,tol=tol,damping_trigger_iter=damping_trigger_iter,
       damping_active=damping_active,trace=do.call(rbind,rows[seq_len(iter)]))
}

## Exactly n_steps raw target omega sweeps; alpha/gamma stay frozen.
## Endpoint convergence is diagnostic only and never gates the next stage.
trace_target_only_fixed_steps <- function(dat, omega, S0bar_inv,
                                          rho_z, A_y=1, B_y=1,
                                          n_steps=5L, tol=1e-6,
                                          verbose=FALSE) {
  stopifnot(length(omega)==dat$p, length(n_steps)==1L,is.finite(n_steps),
            n_steps>=1,n_steps==as.integer(n_steps),
            length(tol)==1L,is.finite(tol),tol>0,
            is.finite(rho_z),rho_z>0,rho_z<1)
  n_steps <- as.integer(n_steps)
  initial_omega <- omega
  rows <- vector("list",n_steps)
  target_map <- function(w) .update_omega_all_woodbury(
    omega=w,gamma=rep(0,dat$p),bbar=rep(0,dat$p),dat=dat,
    S0bar_inv=S0bar_inv,A_y=A_y,B_y=B_y,E_logit_rho_z=qlogis(rho_z))
  for(iter in seq_len(n_steps)) {
    proposed <- target_map(omega)
    if(!all(is.finite(proposed)))stop("Nonfinite fixed-step target update.")
    delta <- max(abs(proposed-omega))
    omega <- proposed
    rows[[iter]] <- data.frame(iter=iter,max_omega_change=delta,
      raw_fixed_point_residual_before=delta,endpoint_fixed_point_residual=NA_real_,
      damping_active=FALSE,damping_factor=1,
      size05=sum(omega>.5),sum_omega=sum(omega),gamma_max=0,alpha_updated=FALSE)
    if(verbose)cat(sprintf("[target fixed] sweep=%d/%d | change=%.4e | size05=%d | gamma=0\n",iter,n_steps,delta,sum(omega>.5)))
  }
  # One extra evaluation measures F(omega_N)-omega_N; its proposal is discarded.
  endpoint_proposed <- target_map(omega)
  if(!all(is.finite(endpoint_proposed)))stop("Nonfinite fixed-step target endpoint residual.")
  residual <- max(abs(endpoint_proposed-omega))
  rows[[n_steps]]$endpoint_fixed_point_residual <- residual
  list(omega=omega,initial_omega=initial_omega,converged=residual<tol,
       phase_complete=TRUE,mode="fixed_steps",n_iter=n_steps,n_steps=n_steps,
       map_evaluations=n_steps+1L,max_change=residual,tol=tol,
       damping_trigger_iter=NA_integer_,damping_active=FALSE,
       trace=do.call(rbind,rows))
}
