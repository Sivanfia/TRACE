# Base-R scientific figures; uncertainty is Monte Carlo SE within each fixed cell.
sim_make_plots <- function(summary, result_root, setting) {
  plot_dir <- file.path(result_root, "figures")
  dir.create(plot_dir, showWarnings = FALSE)
  if (!nrow(summary)) {
    writeLines("No completed replicates; no figures available.", file.path(plot_dir, "README.txt"))
    return(invisible(NULL))
  }
  metrics <- c(f1 = "Variable selection", precision = "Support selection: precision",
    recall = "Support selection: recall", support_exact = "Exact support recovery probability",
    model_size = "Selected support size", borrowing_f1 = "Borrowing set: F1",
    borrowing_precision = "Borrowing set: precision", borrowing_recall = "Borrowing set: recall",
    borrowing_exact = "Exact borrowing-set recovery probability",
    alpha_l1_error = "Source weights: L1 error", pred_mse = "Test-response prediction MSE",
    converged = "Algorithm convergence fraction")
  palette <- c("TRACE-EM" = "#0072B2", "Oracle TRACE" = "#009E73", "ACVB" = "#D55E00",
    "Trans-Lasso" = "#CC79A7", "TRADER" = "#666666", "CONCERT" = "#E69F00")
  symbols <- setNames(c(16, 17, 15, 18, 8, 4), names(palette))
  summary <- summary[summary$metric %in% names(metrics), , drop = FALSE]
  write.csv(summary, file.path(plot_dir, "plot_data.csv"), row.names = FALSE)
  facets <- c("nzc_grid", "rho_grid", "snr_grid", intersect(c("shared_spur_count_grid","sigma_fixed_grid"), names(summary)))
  summary$facet <- do.call(paste, c(summary[facets], sep = "|"))
  index <- list()
  for (metric in names(metrics)) {
    d <- summary[summary$metric == metric & is.finite(summary$mean), , drop = FALSE]
    if (!nrow(d)) next
    d <- d[do.call(order, d[facets]), , drop = FALSE]
    panels <- unique(d$facet); pages <- split(panels, ceiling(seq_along(panels) / 6L))
    methods <- names(palette)[names(palette) %in% unique(d$method)]
    se <- ifelse(is.finite(d$se), d$se, 0)
    proportions <- metric %in% c("f1", "precision", "recall", "support_exact", "borrowing_f1",
      "borrowing_precision", "borrowing_recall", "borrowing_exact", "converged")
    ylim <- range(c(0, if (proportions) 1 else numeric(), d$mean - se, d$mean + se))
    span <- diff(ylim); if (span == 0) span <- 1
    ylim <- ylim + c(-.04, .06) * span
    draw_panel <- function(facet) {
        z <- d[d$facet == facet, , drop = FALSE]; g <- z[1, ]
        if (identical(as.character(setting), "04")) {
          main <- sprintf(
            "s=%d   rho=%g",
            g$nzc_grid, g$rho_grid
          )
        } else {
          main <- sprintf(
            "s=%d   rho=%g   SNR=%.3g",
            g$nzc_grid, g$rho_grid, g$snr_grid
          )
        }
        if ("sigma_fixed_grid" %in% names(g)) main <- paste0(main, "; sigma=", g$sigma_fixed_grid)
        if ("shared_spur_count_grid" %in% names(g)) main <- paste0(main, "\nshared nulls=", g$shared_spur_count_grid)
        xr <- range(z$common_bad_frac_grid)
        if (diff(xr) == 0) xr <- xr + c(-.05, .05)
        plot(NA, xlim = xr + c(-.025, .025), ylim = ylim, xlab = "Contamination fraction q",
          ylab = if(metric=="f1") expression(F[1]) else unname(metrics[metric]), main = main, cex.main = .91, xaxt = "n")
        axis(1, at = sort(unique(z$common_bad_frac_grid)))
        grid(col = "#e8e8e8")
        if (metric == "model_size") abline(h = g$nzc_grid, lty = 3, col = "#999999")
        for (method in methods) {
          x <- z[z$method == method, , drop = FALSE]
          x <- x[order(x$common_bad_frac_grid), , drop = FALSE]
          if (!nrow(x)) next
          ok <- is.finite(x$se) & x$se > 1e-12*pmax(1,abs(x$mean))
          if(any(ok)) {
            xx<-x$common_bad_frac_grid[ok];lo<-x$mean[ok]-x$se[ok];hi<-x$mean[ok]+x$se[ok]
            cap<-(diff(xr)+.05)*.012;ecol<-adjustcolor(palette[method],.6)
            segments(xx,lo,xx,hi,col=ecol,lwd=.8)
            segments(xx-cap,lo,xx+cap,lo,col=ecol,lwd=.8)
            segments(xx-cap,hi,xx+cap,hi,col=ecol,lwd=.8)
          }
          lines(x$common_bad_frac_grid, x$mean, col = palette[method], lwd = 1.6,
            type = "b", pch = symbols[method], cex = .7, lty = match(method, names(palette)))
        }
    }
    draw_page <- function(page_id) {
      panel_ids <- pages[[page_id]]
      ncols <- min(3L, length(panel_ids)); nrows <- ceiling(length(panel_ids) / ncols)
      par(mfrow = c(nrows, ncols), mar = c(3.4, 3.8, 2.8, .8), oma = c(5.2, .3, 3.1, .3),
        mgp = c(2.1, .65, 0), cex = .85)
      for (facet in panel_ids) draw_panel(facet)
      for (unused in seq_len(nrows * ncols - length(panel_ids))) plot.new()
      mtext(sprintf("Setting %s | %s | page %d/%d", setting, metrics[metric], page_id, length(pages)),
        side = 3, outer = TRUE, line = 1.1, cex = 1.15, font = 2)
      nvals <- d$n_eff[d$facet %in% panel_ids]
      count <- if (min(nvals) == max(nvals)) as.character(min(nvals)) else paste(range(nvals), collapse = "-")
      note <- sprintf("Mean +/- 1 Monte Carlo SE; n=%s completed reps per cell (requested %d). SE undefined for n=1.", count, max(d$n_requested))
      mtext(note, side = 1, outer = TRUE, line = 2.5, cex = .83)
      note2 <- if (any(d$alpha_identification_boundary_grid %in% TRUE)) "q=1: C0 is empty; source-weight identification assumptions do not apply."
        else "Each panel holds s, rho and SNR fixed; random seeds are paired across q."
      mtext(note2, side = 1, outer = TRUE, line = 3.65, cex = .8)
      par(fig = c(0, 1, .065, .11), new = TRUE, mar = rep(0, 4), oma = rep(0, 4))
      plot.new()
      legend("center", legend = methods, col = palette[methods],
        pch = symbols[methods], lty = match(methods, names(palette)), lwd = 1.6,
        ncol = length(methods), bty = "n", cex = .92, xpd = NA)
    }
    # Every panel is also a standalone publication-friendly vector PDF / 300dpi PNG.
    single_dir <- file.path(plot_dir,"individual",metric)
    dir.create(single_dir,recursive=TRUE,showWarnings=FALSE)
    token <- function(x) gsub("[^A-Za-z0-9_-]","p",format(x,digits=6,trim=TRUE,scientific=FALSE))
    for(panel_no in seq_along(panels)) {
      facet <- panels[panel_no]; z <- d[d$facet==facet,,drop=FALSE]; g <- z[1,]
      stem <- paste0("setting",setting,"_",metric,"_panel",sprintf("%03d",panel_no),
        "_s",g$nzc_grid,"_rho",token(g$rho_grid),"_snr",token(g$snr_grid))
      if("sigma_fixed_grid"%in%names(g))stem<-paste0(stem,"_sigma",token(g$sigma_fixed_grid))
      if("shared_spur_count_grid"%in%names(g))stem<-paste0(stem,"_spur",g$shared_spur_count_grid)
      draw_single <- function(){
        par(mfrow=c(1,1),mar=c(4.1,4.6,3.1,1),oma=c(4.4,0,1.6,0),mgp=c(2.6,.75,0),cex=.95)
        draw_panel(facet)
        mtext(paste0("Setting ",setting),side=3,outer=TRUE,line=.2,cex=1,font=2)
        nvals<-z$n_eff;count<-if(min(nvals)==max(nvals))as.character(min(nvals))else paste(range(nvals),collapse="-")
        mtext(paste0("Mean +/- 1 Monte Carlo SE; ",count," completed reps/cell. SE undefined for n=1."),
          side=1,outer=TRUE,line=2.2,cex=.72)
        mtext("q is requested contamination fraction; actual bad count = round(s*q).",
          side=1,outer=TRUE,line=3.25,cex=.72)
        par(fig=c(0,1,.102,.165),new=TRUE,mar=rep(0,4),oma=rep(0,4));plot.new()
        legend("center",methods,col=palette[methods],pch=symbols[methods],
          lty=match(methods,names(palette)),lwd=1.6,ncol=min(3L,length(methods)),bty="n",cex=.82)
      }
      local({pdf(file.path(single_dir,paste0(stem,".pdf")),width=7,height=6,useDingbats=FALSE)
        on.exit(dev.off());draw_single()})
      local({type<-if(identical(Sys.info()[["sysname"]],"Darwin"))"quartz"else"cairo"
        png(file.path(single_dir,paste0(stem,".png")),width=2100,height=1800,res=300,type=type)
        on.exit(dev.off());draw_single()})
      index[[length(index)+1L]]<-data.frame(layout="individual",setting=setting,metric=metric,
        nzc=g$nzc_grid,rho=g$rho_grid,snr=g$snr_grid,
        sigma_fixed=if("sigma_fixed_grid"%in%names(g))g$sigma_fixed_grid else NA_real_,
        shared_spur_count=if("shared_spur_count_grid"%in%names(g))g$shared_spur_count_grid else NA_real_,
        pdf=file.path("individual",metric,paste0(stem,".pdf")),page=1L,
        png=file.path("individual",metric,paste0(stem,".png")))
    }
    pdf_name <- paste0(metric, ".pdf")
    local({
      pdf(file.path(plot_dir, pdf_name), width = 12, height = 8, onefile = TRUE, useDingbats = FALSE)
      on.exit(dev.off())
      for (page_id in seq_along(pages)) draw_page(page_id)
    })
    for (page_id in seq_along(pages)) {
      png_name <- sprintf("%s_page%02d.png", metric, page_id)
      local({
        # Quartz on macOS; Cairo on a headless Linux server. Both are standard R devices.
        type <- if (identical(Sys.info()[["sysname"]], "Darwin")) "quartz" else "cairo"
        png(file.path(plot_dir, png_name), width = 1800, height = 1200, res = 150, type = type)
        on.exit(dev.off())
        draw_page(page_id)
      })
      index[[length(index) + 1L]] <- data.frame(layout="combined",setting=setting,metric=metric,nzc=NA,rho=NA,snr=NA,
        sigma_fixed=NA,shared_spur_count=NA,pdf=pdf_name,page=page_id,png=png_name)
    }
  }
  figure_index <- do.call(rbind,index)
  write.csv(figure_index,file.path(plot_dir,"figure_index.csv"),row.names=FALSE)
  write.csv(figure_index[figure_index$layout=="individual",],file.path(plot_dir,"individual_figure_index.csv"),row.names=FALSE)
  writeLines(c("Every panel is saved under individual/<metric>/ as vector PDF (7 x 6 inches) and PNG (2100 x 1800, 300 dpi).",
    "individual_figure_index.csv identifies each panel by setting, s, rho, SNR, sigma and shared-null count.",
    "Combined overview PDF/PNG are also retained. All plots use figures/plot_data.csv.",
    "Error bars are mean +/- one Monte Carlo standard error, not confidence intervals.",
    "Statistics are computed within each fixed scenario cell, never pooled over rho/SNR/q.",
    "Unsupported metrics stay NA. Current selection convention sets empty-denominator precision/recall/F1 to zero.",
    "See n_eff, n_completed and n_requested in the CSV files before interpreting partial results.",
    "Exact borrowing recovery compares C_hat = S_hat intersect Ccred_hat to true C0.",
    "Prediction MSE compares predictions to independent observed test responses, including test noise."),
    file.path(plot_dir, "README.txt"))
  invisible(plot_dir)
}
