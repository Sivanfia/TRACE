## One publication-ready figure per metric and paper setting.
sim_make_plots <- function(summary, result_root, setting) {
  plot_dir <- file.path(result_root, "figures", "individual")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(summary)) {
    writeLines("No completed replicates; no figures available.",
               file.path(plot_dir, "README.txt"))
    return(invisible(NULL))
  }
  labels <- c(
    f1 = "Support F1", precision = "Support precision", recall = "Support recall",
    support_exact = "Exact support recovery probability",
    pred_mse = "Test-response prediction MSE", model_size = "Selected support size",
    alpha_informative_mass = "Estimated mass on informative sources",
    alpha_uninformative_mass = "Estimated mass on uninformative sources",
    time_sec = "Runtime (seconds)", converged = "Convergence fraction"
  )
  palette <- c(
    "TRACE-EM" = "#0072B2", "Oracle TRACE" = "#009E73", "ACVB" = "#D55E00",
    "Trans-Lasso" = "#CC79A7", "TRADER" = "#555555"
  )
  symbols <- setNames(c(16, 17, 15, 18, 8), names(palette))
  ltypes <- setNames(seq_along(palette), names(palette))
  d0 <- summary[summary$metric %in% names(labels), , drop = FALSE]
  sim_csv(d0, file.path(dirname(plot_dir), "plot_data.csv"))
  index <- list()

  draw_plot <- function(d, metric, paper_label) {
    methods <- names(palette)[names(palette) %in% unique(d$method)]
    errs <- ifelse(is.finite(d$se), d$se, 0)
    prop_metric <- metric %in% c(
      "f1", "precision", "recall", "support_exact", "converged",
      "alpha_informative_mass", "alpha_uninformative_mass"
    )
    yr <- range(c(d$mean - errs, d$mean + errs,
                  if (prop_metric) c(0, 1) else numeric()), finite = TRUE)
    if (!all(is.finite(yr))) yr <- c(0, 1)
    if (diff(yr) == 0) yr <- yr + c(-0.05, 0.05)
    yr <- yr + c(-0.04, 0.07) * diff(yr)
    first <- d[1, ]
    plot(NA, xlim = c(1, 10), ylim = yr, xaxs = "i",
      xlab = expression(K[a]), ylab = unname(labels[metric]), xaxt = "n",
      main = sprintf("Setting %s: n=%d, p=%d, s=%d, h=%g, SNR=%.3f",
        paper_label, first$n_grid, first$p_grid, first$nzc_grid, first$h_grid,
        first$snr_grid))
    axis(1, at = 1:10); grid(col = "#E8E8E8")
    if (metric == "model_size") abline(h = first$nzc_grid, lty = 3, col = "#888888")
    for (method in methods) {
      z <- d[d$method == method, , drop = FALSE]
      z <- z[order(z$Ka_grid), , drop = FALSE]
      ok <- is.finite(z$mean)
      if (!any(ok)) next
      z <- z[ok, , drop = FALSE]
      se <- ifelse(is.finite(z$se), z$se, 0)
      has_se <- se > 0
      if (any(has_se)) {
        x <- z$Ka_grid[has_se]; lo <- z$mean[has_se] - se[has_se]
        hi <- z$mean[has_se] + se[has_se]
        ecol <- grDevices::adjustcolor(palette[method], alpha.f = 0.55)
        segments(x, lo, x, hi, col = ecol, lwd = 0.8)
        segments(x - 0.08, lo, x + 0.08, lo, col = ecol, lwd = 0.8)
        segments(x - 0.08, hi, x + 0.08, hi, col = ecol, lwd = 0.8)
      }
      lines(z$Ka_grid, z$mean, type = "b", col = palette[method],
        pch = symbols[method], lty = ltypes[method], lwd = 1.7, cex = 0.8)
    }
    legend("topright", legend = methods, col = palette[methods], pch = symbols[methods],
      lty = ltypes[methods], lwd = 1.7, bty = "n", cex = 0.82)
    mtext("Mean +/- 1 Monte Carlo SE; spike-and-slab source summaries are shared by TRACE and TRADER.",
          side = 1, outer = TRUE, line = 2.1, cex = 0.72)
  }

  for (paper_label in c("I", "III")) {
    for (metric in names(labels)) {
      d <- d0[d0$paper_setting_grid == paper_label & d0$metric == metric &
                is.finite(d0$mean), , drop = FALSE]
      if (!nrow(d)) next
      stem <- sprintf("setting%s_%s_%s", setting, paper_label, metric)
      pdf_file <- file.path(plot_dir, paste0(stem, ".pdf"))
      png_file <- file.path(plot_dir, paste0(stem, ".png"))
      grDevices::pdf(pdf_file, width = 7, height = 5.8, useDingbats = FALSE)
      par(mar = c(4.2, 4.6, 3.2, 1), oma = c(3.2, 0, 0, 0), mgp = c(2.7, 0.75, 0))
      draw_plot(d, metric, paper_label)
      grDevices::dev.off()
      png_type <- if (identical(Sys.info()[["sysname"]], "Darwin")) "quartz" else "cairo"
      grDevices::png(png_file, width = 2100, height = 1740, res = 300, type = png_type)
      par(mar = c(4.2, 4.6, 3.2, 1), oma = c(3.2, 0, 0, 0), mgp = c(2.7, 0.75, 0))
      draw_plot(d, metric, paper_label)
      grDevices::dev.off()
      index[[length(index) + 1L]] <- data.frame(
        setting = setting, paper_setting = paper_label, metric = metric,
        pdf = file.path("individual", basename(pdf_file)),
        png = file.path("individual", basename(png_file)), stringsAsFactors = FALSE
      )
    }
  }
  idx <- if (length(index)) do.call(rbind, index) else data.frame()
  sim_csv(idx, file.path(dirname(plot_dir), "individual_figure_index.csv"))
  writeLines(c(
    "Each metric and paper setting is stored as a separate vector PDF and 300-dpi PNG.",
    "The x-axis is the number Ka of informative sources among K=10.",
    "Error bars are mean +/- one Monte Carlo standard error; they are absent for rep=1.",
    "Coefficient MSE is intentionally not computed or plotted."
  ), file.path(dirname(plot_dir), "README.txt"))
  invisible(plot_dir)
}
