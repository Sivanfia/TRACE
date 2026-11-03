# CONCERT upstream snapshot

This directory contains the linear-model files from the authors' public
repository:

- Repository: https://github.com/RuqianZhang/CONCERT
- Commit: `fea55703c39655a196d3aa5296d1851f0fd53547`
- Retrieved: 2026-09-16

`Linear/Concert.cpp` is compiled and called without algorithmic edits. The
project-level adapter is `../../concert_simulation.R`; it supplies absolute paths,
reproduces the authors' ridge-agreement rule for estimating `q_k`, validates dimensions
and finite output, and maps the posterior mean to the common evaluation
interface.

No license file was present at this upstream commit. Preserve the authors'
attribution and check reuse terms before redistributing these third-party
files outside research reproduction.
