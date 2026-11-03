# Redraw figures without fitting or compiling C++.
a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
p <- if (length(a)) sub("^--file=", "", a[1]) else sys.frame(1)$ofile
if (!file.exists(p)) p <- gsub("~+~", " ", p, fixed = TRUE)
project <- dirname(normalizePath(p))
source(file.path(project, "..", "_common", "plot_results.R"))
root <- path.expand(Sys.getenv("SIM_RESULT_ROOT", file.path(dirname(project), paste0("results_same_source_", Sys.getenv("SIM_SOURCE_N", "500"), "_spikeslab"), basename(project))))
s <- read.csv(file.path(root, "summary_long_all_scenarios.csv"), stringsAsFactors = FALSE)
sim_make_plots(s, root, substr(basename(project), 1, 2))
