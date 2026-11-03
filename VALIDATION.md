# Source distribution checks

Prepared on 3 November 2026 from the complete `simulation_and_real_data_analysis` reproducibility package assembled on 28 October 2026.

## Source preservation

- All 144 original R, C/C++, and shell source files selected for this distribution are byte-identical to the complete package. Their reference MD5 hashes are in `provenance/release_source_md5.csv`.
- Original code comments, method attribution, upstream notes, and environment records are retained.
- The additions are documentation, `.gitignore`, source-distribution manifests, and the data-restoration helper with its test. No statistical algorithm, parameter grid, fitting rule, or plotting calculation was changed.
- External input data, saved fits, results, and manuscript snapshots are omitted from the source archive and remain in the complete Release.

## Checks run

Checks used R 4.6.0 on macOS. A temporary copy was tested outside the source distribution directory.

1. **Restoration helper:** passed fixtures for copying runtime inputs, repeat restoration, rejecting conflicting local data before copying, preserving source files, and rejecting a mismatched source version.
2. **Actual Release restoration:** restored 6,547 runtime files from the complete package. This supplied the inputs and saved records required by the existing reconstruction and validation workflows.
3. **Original validation suite:** `bash run_analysis.sh validate` passed. It checked R syntax, four unchanged grids, simulation and real-data record counts, baseline merging, 108 simulation residual-identity checks, 24 real-data residual-identity checks, six tissue input-reconstruction checks, saved-data fingerprints, and checkpoint safeguards.
4. **Figure and table reconstruction:** `bash run_analysis.sh figures` passed. All 41,100 simulation records and 1,500 real-data records were aggregated. The four simulation summaries and the real-data numerical tables matched the archived paper reference. Main and supplementary figures were rebuilt.
5. **Packaging:** shell scripts passed syntax checks; source hashes, archive contents, file sizes, upload batches, and Git exclusion rules were checked.

The full 50-replication fitting experiments were not rerun for this packaging task. The broader optional `FULL_VALIDATION=1` small-fit suite was not rerun; historical complete-package validation is preserved in `provenance/full_package_validation.md`.

## Archive integrity

`SHA256SUMS` covers the files shipped in this source distribution, excluding the checksum file itself. Verify after extraction with:

```bash
# Linux
sha256sum -c SHA256SUMS

# macOS
shasum -a 256 -c SHA256SUMS
```

The complete Release archive has its own separate checksum list. Do not replace it with this source-only list.
