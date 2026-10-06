# Artandi Lab — passage comparisons and measurement sensitivity

[Open the demo](https://liuzhitong330.github.io/artandi-lab-demo/)

An exploratory, AI-assisted project by Cathy Liu using public Artandi-team data. It supports experimental interpretation, review and record planning; it is not a clinically or experimentally validated assay and does not establish hands-on experience.

## What is here

- A row-preserving analysis of 2,610 deposited TERT-knockout telomere measurements across four passages. Compare mean, median and short-tail fraction, with explicit read-selection stress tests.
- An actual execution of Telometer 2.0.3 on the authors’ public example BAM, comparing three explicit gap settings and exporting read-level review lists.
- Downloadable comparison records and a blank prospective validation log. No entries imply that a proposed experiment was performed.

The two inputs are separate; the example BAM is not mapped to the passage samples. Read counts are not independent biological replication. Keep the stated provenance, processing and limitations with any reused output.

## Reproduce and test

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the publisher-data build and [TELOMETER_AUDIT.md](TELOMETER_AUDIT.md) for the pinned package/BAM reproduction.

```sh
python3 -m unittest discover -s tests -v
python3 scripts/test_telometer_audit.py
node tests/test_app.cjs
python3 -m http.server 8000
```

The first test command needs the dependencies in `requirements-analysis.txt`; the audit test uses the Python standard library. Then open `http://localhost:8000/`. Browser downloads are generated locally; the website collects no entered data and uses no third-party analytics.

## Sources and reuse

- Sanchez et al. (2024), *Nature Communications*, [doi:10.1038/s41467-024-49007-4](https://doi.org/10.1038/s41467-024-49007-4). Published source-data subset and its transformations retain [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) attribution.
- [santiago-es/Telometer](https://github.com/santiago-es/Telometer), pinned commit `f17a7e43e6d6e4313354180fa42b761b581b1b14`; MIT notice in `licenses/Telometer-MIT.txt`.
- New demo code is licensed under MIT; see `LICENSE`. Upstream data/code retain their own licenses. Raw sequencing reads, private application files, CVs and reference contacts are not included.

Cathy Liu · cathyliu014@gmail.com · 415-216-3799
