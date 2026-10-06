# Telometer example-data audit

This is a separate reproducibility check of the Artandi team's public example
alignment file. It is not a reanalysis of the same samples as the TERT-knockout
passage series. It does not demonstrate experience preparing sequencing libraries.

## Provenance

- Primary paper: Sanchez et al., *Nature Communications* 15, 5148 (2024),
  [doi:10.1038/s41467-024-49007-4](https://www.nature.com/articles/s41467-024-49007-4).
  Its code-availability statement identifies the repository below.
- Repository: [santiago-es/Telometer](https://github.com/santiago-es/Telometer),
  commit `f17a7e43e6d6e4313354180fa42b761b581b1b14` (18 September 2026).
- Input: `minimal_tels.bam`, SHA256
  `1b17bcdc765ae62b07a94cf3b2303666f0cf2892d83f2f72dd28dbe508b76af6`.
  The repository describes this as a subsampled telomere-capture experiment.
- Code: Telometer 2.0.3, the PyPI wheel linked through the repository's installation
  instructions, SHA256 `1314094e09a8cb494521b60cf70caef2d8093638b0d3964be449558fbfc35348`.
  The Python package is newer than the publication. This is therefore a pinned
  current-version test, **not an exact reproduction of the 2024 analysis**.
- The repository and code package use the MIT license. The upstream notice is
  retained in `licenses/Telometer-MIT.txt`. No sequencing reads are redistributed
  with this website; it contains only derived measurements and public read IDs.

## Question

Which example reads should receive manual review when telomere-boundary gap
tolerance changes? This helps distinguish a stable distribution-level result
from a small number of parameter-sensitive measurements before planning an
orthogonal assay. It does not identify the biologically correct gap tolerance.

## Reproduce

Use Python 3.12 and an isolated environment. The checked run used macOS arm64.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r scripts/requirements-telometer.txt
.venv/bin/python scripts/reproduce_telometer_audit.py --cache-dir ../artandi-analysis-cache
```

The script downloads and SHA256-checks the pinned BAM and package wheel, checks
the installed package version, runs three explicit configurations, and writes
`data/telometer-audit.json`. Keep the cache and environment outside the public
repository. The script does not install or execute the downloaded wheel itself;
installation is the explicit preceding command. `--reuse-tsv` is only for results
already generated from exactly these inputs and settings.

Common settings: minimum read length 1,000 bp, two worker processes, one GB batch
target. Gap tolerances are 20, 100 and 250 bp. Version 2.0.3's CLI default is 20 bp;
the pinned README separately gives 100 bp in an example and 250 bp in its options
table. Explicit settings remove that ambiguity. No optional variant motif or
enrichment calculation is requested.

Rows are joined by original read ID. Each setting retains its own selected
chromosome, arm, direction, MAPQ and read length because Telometer's winning
alignment can change as its length measurement changes. An absent output is
stored as `null`, not zero. The included `read-001` labels are stable ordinal
display labels assigned after sorting public read IDs, not biological sample IDs.

## Limits

The example contains 1,640 alignment records from 812 distinct read IDs. Alignment
records are not independent molecules. A Telometer output retains one selected
alignment per read ID, and that does not validate chromosome assignment or prove
an intact telomere. The BAM's sampling, biological-replicate identity and full
sequencing yield are not established by this test. No enrichment, disease status,
cell-state or between-culture biological uncertainty should be inferred from it.

A sensitivity flag means inspect the underlying read and alignment. It does not
mean discard it automatically. A stable median can coexist with unstable extreme
measurements. This audit is a reusable QC prototype, not a validated lab SOP.
