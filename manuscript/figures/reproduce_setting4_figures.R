# Rebuild Setting 4 supplementary diagnostics; keep the common main-figure style.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
if (!file.exists(entry)) entry <- gsub("~+~", " ", entry, fixed=TRUE)
figure_dir <- dirname(normalizePath(entry))
project <- dirname(figure_dir)
latest <- read.csv(file.path(project,"data","setting04_summary.csv"),stringsAsFactors=FALSE)
methods <- c("TRACE-EM","Oracle TRACE","ACVB","Trans-Lasso","TRADER","CONCERT")
colors <- setNames(c("#0072B2","#009E73","#D55E00","#CC79A7","#666666","#E69F00"),methods)
symbols <- setNames(c(16,17,15,18,8,4),methods)
stopifnot(all(latest$n_completed==50),all(latest$n_requested==50),
          all(latest$rho_grid==.3),all(latest$sigma_fixed_grid==3),
          length(unique(latest$scenario_index))==21)

error_lines <- function(x,y,se,col,pch,lty=1,lwd=1.1) {
  good <- is.finite(x)&is.finite(y)&is.finite(se)
  x<-x[good];y<-y[good];se<-se[good]
  order <- order(x);x<-x[order];y<-y[order];se<-se[order]
  cap <- .005
  segments(x,y-se,x,y+se,col=adjustcolor(col,.65),lwd=.6)
  segments(x-cap,y-se,x+cap,y-se,col=adjustcolor(col,.65),lwd=.6)
  segments(x-cap,y+se,x+cap,y+se,col=adjustcolor(col,.65),lwd=.6)
  lines(x,y,type="b",col=col,pch=pch,lty=lty,lwd=lwd,cex=.53)
}

f1_plot <- function(data,sigma,ticks,filename) {
  d <- data[data$metric=="f1"&data$sigma_fixed_grid==sigma,]
  pdf(file.path(figure_dir,filename),width=6.5,height=2.3,pointsize=9,
      family="Helvetica",useDingbats=FALSE)
  layout(rbind(1:3,rep(4,3)),heights=c(1,.25))
  for (i in 1:3) {
    s<-c(10,15,30)[i];z<-d[d$nzc_grid==s,]
    stopifnot(nrow(z)==length(methods)*length(ticks),all(z$n_eff==50))
    par(cex=1,mar=c(2.65,if(i==1)3.3 else .85,1.5,.25),
        mgp=c(1.65,.45,0),tcl=-.2,las=1,cex.axis=.82,cex.lab=.9,
        xaxs="i",yaxs="i")
    pad<-diff(range(ticks))*.055
    plot(NA,xlim=range(ticks)+c(-pad,pad),ylim=c(-.025,1.055),
         xaxt="n",yaxt="n",xlab="Contamination fraction q",
         ylab=if(i==1)expression(F[1]) else "",bty="l")
    abline(h=seq(0,1,.2),col="#E8E8E8",lwd=.55)
    axis(1,at=ticks,labels=sprintf("%.1f",ticks))
    if(i==1)axis(2,at=seq(0,1,.2),labels=sprintf("%.1f",seq(0,1,.2)))
    title(main=bquote(.(paste0("(",letters[i],") "))*italic(s)==.(s)),
          cex.main=1.04,font.main=1,line=.6)
    for(method in rev(methods)) {
      a<-z[z$method==method,]
      stopifnot(all(is.finite(a$mean)),all(is.finite(a$se)))
      error_lines(a$common_bad_frac_grid,a$mean,a$se,colors[method],symbols[method],match(method,methods))
    }
  }
  par(mar=c(0,0,0,0));plot.new()
  legend("center",methods,col=colors,pch=symbols,lty=seq_along(methods),
         lwd=1.1,pt.cex=.65,ncol=3,bty="n",cex=.95,x.intersp=.8,y.intersp=1.05)
  dev.off()
}

diagnostic_plot <- function() {
  d<-latest[latest$method=="TRACE-EM",]
  pdf(file.path(figure_dir,"setting04_transfer_diagnostics.pdf"),
      width=6.5,height=3.7,pointsize=9,family="Helvetica",useDingbats=FALSE)
  layout(rbind(1:3,4:6,rep(7,3)),heights=c(1,1,.23))
  for(row in 1:2)for(i in 1:3) {
    s<-c(10,15,30)[i]
    par(cex=1,mar=c(3.0,if(i==1)4.2 else .85,1.5,.25),
        mgp=c(1.65,.45,0),tcl=-.2,las=1,cex.axis=.82,cex.lab=.9,
        xaxs="i",yaxs="i")
    plot(NA,xlim=c(-.025,.625),ylim=if(row==1)c(-.004,.19) else c(-.025,1.055),
         xaxt="n",yaxt="n",xlab="Contamination fraction q",bty="l",
         ylab="")
    if(i==1)mtext(if(row==1)"Weight L1 error" else "Mean credibility",
                  side=2,line=2.8,cex=.9,las=0)
    yticks<-if(row==1)seq(0,.18,.06)else seq(0,1,.2)
    abline(h=yticks,col="#E8E8E8",lwd=.55)
    axis(1,at=(0:6)/10,labels=sprintf("%.1f",(0:6)/10))
    if(i==1)axis(2,at=yticks,labels=sprintf(if(row==1)"%.2f" else "%.1f",yticks))
    title(main=bquote(.(paste0("(",letters[(row-1)*3+i],") "))*italic(s)==.(s)),
          cex.main=1.04,font.main=1,line=.6)
    if(row==1) {
      a<-d[d$nzc_grid==s&d$metric=="alpha_l1_error",]
      stopifnot(nrow(a)==7,all(a$n_eff==50))
      error_lines(a$common_bad_frac_grid,a$mean,a$se,"#0072B2",16)
    }else {
      for(class in c("mean_gamma_transferable_signal","mean_gamma_common_bad")) {
        a<-d[d$nzc_grid==s&d$metric==class,]
        compatible<-class=="mean_gamma_transferable_signal"
        error_lines(a$common_bad_frac_grid,a$mean,a$se,
                    if(compatible)"#009E73"else"#D55E00",
                    if(compatible)17 else 15,if(compatible)1 else 2)
      }
    }
  }
  par(mar=c(0,0,0,0));plot.new()
  legend("center",c("Weight error","Transferable active","Contaminated active"),
         col=c("#0072B2","#009E73","#D55E00"),pch=c(16,17,15),lty=c(1,1,2),
         lwd=1.1,pt.cex=.65,ncol=3,bty="n",cex=.95,x.intersp=.8)
  dev.off()
}

# Main-text Setting 4 uses the common style in reproduce_simulation_figures.R.
# f1_plot(latest,3,(0:6)/10,"setting04_main_f1.pdf")
diagnostic_plot()
cat("Rebuilt Setting 4 supplementary diagnostics.\n")
