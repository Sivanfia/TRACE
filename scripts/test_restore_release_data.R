#!/usr/bin/env Rscript
# Integration tests for restoring external inputs without replacing source code.
entry <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1L])
root <- dirname(dirname(normalizePath(entry)))
test_restore <- function() {
  stopifnot(file.exists(file.path(root, 'scripts/restore_release_data.R')))
  work <- tempfile('TRACE_restore_test_')
  dir.create(work)
  on.exit(unlink(work, recursive=TRUE), add=TRUE)
  release <- file.path(work, 'release')
  checkout <- file.path(work, 'checkout')
  for (d in c(release, checkout)) {
    dir.create(file.path(d, 'scripts'), recursive=TRUE)
    dir.create(file.path(d, 'provenance'))
  }
  writeLines('# original source', file.path(release, 'scripts/run.R'))
  dir.create(file.path(checkout, 'plotting'))
  file.copy(file.path(release, 'scripts/run.R'), file.path(checkout, 'plotting/plot.R'))
  writeLines('# source-layout adapter', file.path(checkout, 'scripts/run.R'))
  write.csv(data.frame(file=c('plotting/plot.R', 'scripts/run.R'),
    release_file=rep('scripts/run.R', 2L),
    md5=unname(tools::md5sum(file.path(checkout, c('plotting/plot.R', 'scripts/run.R')))),
    release_md5=rep(unname(tools::md5sum(file.path(release, 'scripts/run.R'))), 2L)),
    file.path(checkout, 'provenance/release_source_md5.csv'), row.names=FALSE)
  file.copy(file.path(root, 'scripts/restore_release_data.R'),
    file.path(checkout, 'scripts/restore_release_data.R'))
  writeLines('{"name":"simulation_and_real_data_analysis"}', file.path(release, 'MANIFEST.json'))
  paths <- c('simulation/retained_results/01/reference.rds',
    'simulation/published_results/01_TRACE_identifiable_mixture_contamination/comparison/performance_all_scenarios.csv',
    'real_data/data/processed/gtex/tissue.tsv.gz',
    'real_data/retained_results/jobs/a.rds',
    'real_data/published_results/comparison/jobs/a.rds',
    'real_data/published_results/comparison/tables/metrics_per_rep.csv',
    'manuscript/data/setting01_summary.csv', 'manuscript/draft_revised_crossrefs.tex')
  for (rel in paths) {
    dir.create(dirname(file.path(release, rel)), recursive=TRUE, showWarnings=FALSE)
    writeLines(paste('fixture', rel), file.path(release, rel))
  }
  dir.create(file.path(release, 'manuscript/figures'), recursive=TRUE)
  writeLines('# should not be copied', file.path(release, 'manuscript/figures/plot.R'))
  invoke <- function() {
    x <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'),
      c(shQuote(file.path(checkout, 'scripts/restore_release_data.R')), shQuote(release)),
      stdout=TRUE, stderr=TRUE))
    list(status=if(is.null(attr(x, 'status'))) 0L else attr(x, 'status'), output=x)
  }
  x <- invoke()
  if(x$status != 0L) stop(paste(x$output, collapse='\n'))
  stopifnot(all(file.exists(file.path(checkout, paths))),
    !file.exists(file.path(checkout, 'manuscript/figures/plot.R')),
    identical(readLines(file.path(checkout, 'plotting/plot.R')), '# original source'),
    identical(readLines(file.path(checkout, 'scripts/run.R')), '# source-layout adapter'))
  stopifnot(invoke()$status == 0L) # Repeating an identical restoration is safe.
  writeLines('existing local changes', file.path(checkout, paths[1L]))
  x <- invoke()
  stopifnot(x$status != 0L, any(grepl('conflicting', x$output, ignore.case=TRUE)),
    identical(readLines(file.path(checkout, paths[1L])), 'existing local changes'))
  unlink(file.path(checkout, paths))
  writeLines('# a different source version', file.path(release, 'scripts/run.R'))
  x <- invoke()
  stopifnot(x$status != 0L, any(grepl('source', x$output, ignore.case=TRUE)),
    !any(file.exists(file.path(checkout, paths))))
  cat('PASS: restoration, repeat use, conflict rejection, source-version check\n')
}
test_restore()
