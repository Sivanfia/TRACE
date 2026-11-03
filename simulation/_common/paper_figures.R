paper_make_figures <- function(summary,result_root,setting) {
  d<-summary[summary$metric=='f1' & is.finite(summary$mean),,drop=FALSE]
  if(!nrow(d))return(invisible(NULL))
  out<-file.path(result_root,'paper_figures');dir.create(out,showWarnings=FALSE)
  palette<-c('TRACE-EM'='#0072B2','Oracle TRACE'='#009E73','ACVB'='#D55E00',
             'Trans-Lasso'='#CC79A7','TRADER'='#666666','CONCERT'='#E69F00')
  pch<-setNames(c(16,17,15,18,8,4),names(palette))
  methods<-names(palette)[names(palette)%in%d$method]
  draw_panel<-function(z,title,xname) {
    if(!nrow(z)){plot.new();title(main=title,cex.main=.9);text(.5,.5,'Pending',cex=.8,col='grey50');return(invisible(NULL))}
    x<-if(setting=='02')z$Ka_grid else z$common_bad_frac_grid
    xr<-range(x,finite=TRUE)
    plot(xr,c(0,1.02),type='n',xlab=xname,ylab=expression(F[1]),main=title,
         xaxt='n',yaxs='i',cex.main=.9,cex.lab=.9,cex.axis=.8)
    axis(1,at=sort(unique(x)),cex.axis=.8)
    abline(h=seq(0,1,.2),col='grey92',lwd=.6)
    for(method in methods) {
      q<-z[z$method==method,,drop=FALSE]
      if(!nrow(q))next
      xx<-if(setting=='02')q$Ka_grid else q$common_bad_frac_grid
      ord<-order(xx);xx<-xx[ord];q<-q[ord,,drop=FALSE]
      se<-ifelse(is.finite(q$se),q$se,0)
      lines(xx,q$mean,col=palette[method],lwd=1.2,lty=if(method=='Oracle TRACE')2 else 1)
      positive<-se>0
      if(any(positive))arrows(xx[positive],pmax(0,q$mean[positive]-se[positive]),
             xx[positive],pmin(1,q$mean[positive]+se[positive]),angle=90,code=3,
             length=.025,col=palette[method],lwd=.65)
      points(xx,q$mean,col=palette[method],pch=pch[method],cex=.58)
    }
    box(lwd=.65)
  }
  draw<-function() {
    if(setting%in%c('01','03')) {
      snr<-if(setting=='01')4 else 3
      z<-d[abs(d$snr_grid-snr)<1e-7,,drop=FALSE]
      layout(matrix(c(1:6,rep(7,3)),nrow=3,byrow=TRUE),heights=c(1,1,.28))
      par(mar=c(3,3,1.45,.45),mgp=c(1.8,.55,0),tcl=-.25)
      for(rho in c(.3,.6))for(s in c(10,15,30))
        draw_panel(z[z$nzc_grid==s & abs(z$rho_grid-rho)<1e-7,,drop=FALSE],
                   bquote(s==.(s)~','~rho==.(rho)),expression(q))
    } else if(setting=='02') {
      layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(1,.22))
      par(mar=c(3.1,3,1.3,.55),mgp=c(1.8,.55,0),tcl=-.25)
      for(design in c('I','III'))draw_panel(d[d$paper_setting_grid==design,,drop=FALSE],
        paste('Design',design),expression(K[a]))
    } else {
      layout(matrix(c(1:3,rep(4,3)),nrow=2,byrow=TRUE),heights=c(1,.22))
      par(mar=c(3,3,1.5,.5),mgp=c(1.8,.55,0),tcl=-.25)
      for(s in c(10,15,30))draw_panel(d[d$nzc_grid==s,,drop=FALSE],bquote(s==.(s)),expression(q))
    }
    par(mar=c(0,0,0,0));plot.new()
    legend('center',legend=methods,col=palette[methods],pch=pch[methods],
           lty=ifelse(methods=='Oracle TRACE',2,1),lwd=1.2,ncol=3,bty='n',cex=.83,x.intersp=.75)
  }
  width<-7.5;height<-if(setting%in%c('01','03'))4.35 else 2.65
  name<-paste0('setting',setting,'_main_f1')
  grDevices::pdf(file.path(out,paste0(name,'.pdf')),width=width,height=height,useDingbats=FALSE)
  tryCatch(draw(),finally=dev.off())
  grDevices::png(file.path(out,paste0(name,'.png')),width=width,height=height,units='in',res=250,
                type=if(Sys.info()[['sysname']]=='Darwin')'quartz' else 'cairo')
  tryCatch(draw(),finally=dev.off())
  plotted<-d
  if(setting%in%c('01','03')){
    selected_snr<-if(setting=='01')4 else 3
    plotted<-d[abs(d$snr_grid-selected_snr)<1e-7,,drop=FALSE]
  }
  write.csv(plotted,file.path(out,paste0(name,'_plot_data.csv')),row.names=FALSE)
  invisible(NULL)
}
