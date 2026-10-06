# Reproducing the TERT-loss passage analysis

This is a descriptive reanalysis of **processed molecule-level telomere measurements**, not a reanalysis of nanopore raw signal. The question is whether the choice of passage pair, summary statistic or read-retention rule changes the follow-up interpretation.

## Source and attribution

Sanchez et al. (2024), *Nature Communications*, [doi:10.1038/s41467-024-49007-4](https://doi.org/10.1038/s41467-024-49007-4). The data come from `Source Data File.xlsx`, sheet `hesc_tertko`, supplied with the paper under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

[Download the publisher's source-data archive](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-49007-4/MediaObjects/41467_2024_49007_MOESM10_ESM.zip).

- Archive SHA256: `9526cdddcbbf793ba0e17fbc9007c82d2f50073e31b0e946e1de4e97c3e26c9b`
- Workbook SHA256: `26115393784129d05c68d92138c4a417df1b5cee9a308976065e37cdbcff4d24`
- Data rows: physical Excel rows 3–2612, all 2,610 rows retained.
- Original fields: `chromosome`, `telomere_length`, `read_length`, `arm`, `direction`, `day`, `Genotype`.
- All source genotypes are `TERT KO`. The normalized subset retains all variable fields and assigns physical worksheet row as a stable audit key. That key is **not** an original molecule ID.

The source-derived subset and analyses are a transformation by Cathy Liu with AI assistance. Source data retain their CC BY 4.0 attribution. No source article figures or images are redistributed.

## Commands

Use Python 3.10 or newer. From this repository's root:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python scripts/build_tertko.py
.venv/bin/python -m unittest discover -s tests -v
```

With Node.js 20 or newer, the interface's exported arithmetic can also be checked offline:

```sh
node tests/test_app.cjs
```

This checks all 234 passage/threshold/read-length settings exposed by the page, including 108 cross-checks against the separately generated Python sensitivity matrix. It also checks all 12 Telometer review configurations, the difference between zero and missing output, threshold boundaries, omission of reads absent from both compared runs, and CSV escaping. The review-list tests use the separate executed audit file described below; regenerate that audit before rerunning those tests if it is absent.

The build downloads the pinned publisher archive, checks the archive hash, reads only the exact workbook member in memory, checks its hash, and regenerates the three files below. It never executes downloaded code or extracts archive paths onto the filesystem. If the source changes, the build stops rather than silently substituting a different dataset.

For an already-downloaded workbook, the same hash check applies:

```sh
.venv/bin/python scripts/build_tertko.py --workbook '/path/to/Source Data File.xlsx'
```

Outputs:

- `data/tertko.json`: row-level data, source notes, day-level summaries and units in field names.
- `data/tertko.csv`: the identical normalized row-level subset for reuse.
- `data/analysis_checks.json`: source audit, repeated visible-row groups and a 144-case sensitivity matrix.

## Descriptive calculations

Mean is the arithmetic mean. Median is the midpoint for an even-sized ordered sample. Quartiles use linear interpolation at `(n − 1) × p`, the same definition as R type 7 and NumPy's default linear quantile. The short-tail fraction is `count(telomere_bp < threshold) / retained molecules`; equality to the threshold is excluded. Read retention uses `read_bp >= minimum_read_bp`.

The audit evaluates all six chronological passage pairs, thresholds 1,000 / 1,500 / 2,000 / 2,500 / 3,000 / 4,000 bp, and minimum read lengths 0 / 1,000 / 4,000 / 8,000 bp. Zero means all deposited rows. The 1,000-bp test is identical to the unfiltered data because the shortest deposited read is 1,117 bp. This is an implementation check, not evidence that a 1,000-bp cutoff is biologically appropriate. The page can additionally interpolate threshold choices in 250-bp steps directly from the same rows.

These filters are intentionally labeled **sensitivity tests**, not recommended QC. Read and telomere lengths are coupled, so restricting read length can change the sampled telomere distribution. A changed result cannot establish that one filter is more accurate or that a technical artifact caused the original result.

## Source audit and checked findings

| Deposited day | Molecules | Mean, bp | Median, bp | Below 2,000 bp |
| --- | ---: | ---: | ---: | ---: |
| 66 | 377 | 4,296.91 | 3,886 | 32 / 377 (8.49%) |
| 78 | 350 | 3,345.61 | 2,824.5 | 95 / 350 (27.14%) |
| 98 | 968 | 2,617.24 | 2,251 | 409 / 968 (42.25%) |
| 105 | 915 | 2,639.39 | 2,066 | 423 / 915 (46.23%) |

From day 98 to deposited day 105, the mean rises **22.15 bp**, while the median falls **185 bp** and the fraction below 2 kb rises **3.98 percentage points**. This is a descriptive endpoint disagreement, not evidence of telomerase reactivation, an identified mechanism, or a statistically established biological effect.

All deposited lengths are finite and positive. No telomere length exceeds its associated read length. The 25 repeated seven-field tuples are retained because this sheet has no molecule ID that could distinguish repeated molecules from distinct molecules with identical measurements. Repeated visible row groups are listed in the audit output. There is no imputation, deduplication or hidden exclusion.

The workbook and article narrative identify the final sample as **105 days**, but the Figure 2 caption says **108 days**. The demo preserves the deposited value and flags the discrepancy; it does not silently resolve it.

## Interpretation boundaries and next experiment

The rows are molecules nested within four deposited passage distributions. They are **not 2,610 independent biological replicates**. Clone identities and biological replicate identifiers are not available in this sheet. The demo therefore does not calculate biological confidence intervals, hypothesis-test p-values, causal effects, shortening rates per division, or experimental power from the molecule count.

To investigate an endpoint disagreement, prospectively repeat the relevant comparison using independent edited clones and matched controls. Record genotype, passage and population doubling, extraction batch, library batch, capture conditions and read depth, and use an orthogonal telomere-length measurement. These are proposed controls, not experiments already completed. The data cannot determine the appropriate number of independent clones without additional variance and design information.

The separate public example-BAM/Telometer reproduction is not matched to these four passage samples. Do not combine its read counts with this dataset or use its subsetted BAM to estimate enrichment. See the separate tool-reproduction record for its exact version, input and commands.

The demo supports experimental planning, assay interpretation and transparent research records. It does not demonstrate independent hESC culture, genome editing, nanopore library preparation or telomerase assay execution.
