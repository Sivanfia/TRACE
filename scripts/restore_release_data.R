#!/usr/bin/env Rscript
# Restore runtime payload from the extracted complete Release, never source files.
entry <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1L])
if (!file.exists(entry)) entry <- gsub('~+~', ' ', entry, fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
args <- commandArgs(TRUE)
if (length(args) != 1L) stop(
  'Usage: Rscript scripts/restore_release_data.R /path/to/simulation_and_real_data_analysis')
release <- normalizePath(args[1L], mustWork=TRUE)
if (identical(release, root)) stop('Extract the complete Release into a separate directory.')
manifest <- file.path(release, 'MANIFEST.json')
if (!file.exists(manifest) || !any(grepl(
    '"name"[[:space:]]*:[[:space:]]*"simulation_and_real_data_analysis"',
    readLines(manifest, warn=FALSE)))) stop('Not a complete TRACE Release directory.')

# Match the archived computational source before associating it with saved data.
source_manifest <- read.csv(file.path(root, 'provenance/release_source_md5.csv'),
  stringsAsFactors=FALSE)
for (base in c(root, release)) {
  paths <- file.path(base, source_manifest$file)
  if (any(!file.exists(paths))) stop('Missing source files in: ', base)
  hashes <- unname(tools::md5sum(paths))
  if (anyNA(hashes) || any(hashes != source_manifest$md5))
    stop('The source version does not match the archived experiment: ', base)
}
collect <- function(relative) {
  directory <- file.path(release, relative)
  if (!dir.exists(directory)) stop('Required Release directory is missing: ', relative)
  paths <- list.files(directory, recursive=TRUE, full.names=TRUE, all.files=TRUE,
    no..=TRUE, include.dirs=FALSE)
  paths[!dir.exists(paths)]
}
payload <- c(
  collect('simulation/retained_results'),
  grep('/comparison/performance_all_scenarios[.]csv$',
    collect('simulation/published_results'), value=TRUE),
  collect('real_data/data'),
  collect('real_data/retained_results'),
  collect('real_data/published_results/comparison/jobs'),
  file.path(release, 'real_data/published_results/comparison/tables/metrics_per_rep.csv'),
  collect('manuscript'))
payload <- sort(unique(payload))
# The manuscript directory contains plotting sources already supplied in main.
payload <- payload[!grepl('[.](R|r|cpp|c|h|hpp|sh)$', payload)]
if (any(!file.exists(payload))) stop('The Release is missing a required result table.')
relative <- substring(payload, nchar(release) + 2L)
destinations <- file.path(root, relative)
exists <- file.exists(destinations)
if (any(exists)) {
  source_hash <- unname(tools::md5sum(payload[exists]))
  destination_hash <- unname(tools::md5sum(destinations[exists]))
  conflicts <- which(exists)[is.na(destination_hash) | source_hash != destination_hash]
  if (length(conflicts)) stop('Conflicting local data; no files copied. Move or rename: ',
    paste(head(relative[conflicts], 5L), collapse=', '))
}
pending <- which(!exists)
cat('Restoring', length(pending), 'files; keeping', sum(exists), 'identical files.\n')
for (i in pending) {
  dir.create(dirname(destinations[i]), recursive=TRUE, showWarnings=FALSE)
  if (!file.copy(payload[i], destinations[i], overwrite=FALSE, copy.date=TRUE))
    stop('Could not copy: ', relative[i])
}
cat('Ready. Run: bash run_analysis.sh figures\n')
