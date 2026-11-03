# Rebuild every simulation figure from the project-contained summaries.
a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
entry <- sub("^--file=", "", a[1])
if (!file.exists(entry)) entry <- gsub("~+~", " ", entry, fixed = TRUE)
figures <- dirname(normalizePath(entry, mustWork = TRUE))
project <- dirname(figures)
if (basename(figures) == "plotting") project <- file.path(project, "manuscript")
source(file.path(figures, "reproduce_simulation_figures.R"))
source(file.path(figures, "reproduce_complete_grids.R"))
for (id in c("01", "03")) {
  d <- read.csv(file.path(project, "data", paste0("setting", id, "_summary.csv")), stringsAsFactors = FALSE)
  d <- d[d$metric == "f1", , drop = FALSE]
  output <- tempfile(paste0("TRACE_figures_", id, "_"))
  dir.create(output)
  sim_make_plots(d, output, id)
  stopifnot(file.copy(file.path(output, "figures", "f1.pdf"),
                     file.path(figures, paste0("setting", id, "_all_f1.pdf")), overwrite = TRUE))
  unlink(output, recursive = TRUE)
}

source(file.path(figures, "reproduce_setting4_figures.R"))

source(file.path(figures,"reproduce_real_data_figures.R"))
