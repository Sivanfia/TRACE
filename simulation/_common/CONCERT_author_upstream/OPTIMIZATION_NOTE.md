# Local CONCERT performance patch

This directory vendors upstream CONCERT commit
`fea55703c39655a196d3aa5296d1851f0fd53547`.

`Linear/Concert.cpp` has been locally refactored for computational efficiency
without changing the statistical CAVI update equations or their Gauss-Seidel
coordinate order. The original vendored source is retained as
`Linear/Concert.cpp.original`.

Main implementation changes:

- target/source row blocks are created once, not inside coordinate loops;
- column squared norms and source `X^T y` vectors are cached once;
- excluded-coordinate products are evaluated from maintained residuals;
- residuals are refreshed from scratch between major sweeps to limit numerical drift;
- variance updates use matrix-vector quadratic forms instead of repeatedly forming `X^T X`;
- iteration/convergence and phase timing diagnostics are returned.

An independent algebraic equivalence check of the old and refactored update
expressions matched all tested posterior quantities to approximately 1e-15 in
double precision over seven sequential CAVI sweeps. The build environment used
to prepare this patch did not contain R, so final `Rcpp::sourceCpp` compilation
must be performed in the target R installation.
