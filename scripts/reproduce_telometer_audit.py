#!/usr/bin/env python3
"""Run pinned Telometer on its public example BAM, then audit gap sensitivity.

This is a 2026 package-version test, not a reproduction of the 2024 paper's
measurements. The BAM is already subsampled. Never calculate enrichment here.
Only the selected TSV output is read; the underlying input is not modified.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import urllib.request

import numpy as np
import pandas as pd
import pysam

REPO = "https://github.com/santiago-es/Telometer"
COMMIT = "f17a7e43e6d6e4313354180fa42b761b581b1b14"
BAM_SHA = "1b17bcdc765ae62b07a94cf3b2303666f0cf2892d83f2f72dd28dbe508b76af6"
WHEEL_SHA = "1314094e09a8cb494521b60cf70caef2d8093638b0d3964be449558fbfc35348"
WHEEL_URL = "https://files.pythonhosted.org/packages/9d/c7/418213ff0585e168f47839fd511b4a01e7773fa721fcc1df73abf4dc20ea/telometer-2.0.3-py3-none-any.whl"
GAPS = (20, 100, 250)


def checked_download(url, path, expected):
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"SHA256 mismatch for {path.name}: {actual}")


def summary(frame):
    x = frame.telomere_length.to_numpy()
    return {"n": int(len(x)), "min_bp": int(x.min()), "q1_bp": float(np.quantile(x, .25)),
            "median_bp": float(np.median(x)), "mean_bp": float(np.mean(x)),
            "q3_bp": float(np.quantile(x, .75)), "max_bp": int(x.max()),
            "below_2000_n": int((x < 2000).sum())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True,
                        help="Private/local cache; must not be committed to the site")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data" / "telometer-audit.json")
    parser.add_argument("--reuse-tsv", action="store_true",
                        help="Only reuse previously computed files from these exact pinned inputs/configurations")
    args = parser.parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    if importlib.metadata.version("telometer") != "2.0.3":
        raise ValueError("Install the pinned Telometer 2.0.3 package first")
    bam = args.cache_dir / "minimal_tels.bam"
    checked_download(f"https://raw.githubusercontent.com/santiago-es/Telometer/{COMMIT}/minimal_tels.bam", bam, BAM_SHA)
    wheel = args.cache_dir / "telometer-2.0.3-py3-none-any.whl"
    checked_download(WHEEL_URL, wheel, WHEEL_SHA)
    frames = {}
    for gap in GAPS:
        out = args.cache_dir / f"telomeres-gap{gap}.tsv"
        command = [sys.executable, "-m", "telometer.telometer", "-b", str(bam),
                   "-o", str(out), "-m", "1000", "-g", str(gap), "-t", "2", "-l", "1"]
        if not (args.reuse_tsv and out.exists()):
            subprocess.run(command, check=True)
        frame = pd.read_csv(out, sep="\t")
        if frame.read_id.duplicated().any() or len(frame) == 0:
            raise ValueError("Expected a nonempty one-row-per-selected-read output")
        if (frame.telomere_length <= 0).any() or (frame.telomere_length >= frame.read_length).any():
            raise ValueError("Invalid telomere length in selected output")
        frames[gap] = frame.set_index("read_id")

    with pysam.AlignmentFile(bam, "rb") as fh:
        records = list(fh)
    bam_stats = {
        "alignment_records": len(records),
        "distinct_read_ids": len({r.query_name for r in records}),
        "primary_records": sum(not r.is_secondary and not r.is_supplementary for r in records),
        "secondary_records": sum(r.is_secondary for r in records),
        "supplementary_records": sum(r.is_supplementary for r in records),
    }
    read_ids = sorted(set().union(*(set(f.index) for f in frames.values())))
    rows = []
    for idx, read_id in enumerate(read_ids, 1):
        row = {"id": f"read-{idx:03d}", "read_id": read_id}
        for gap in GAPS:
            if read_id not in frames[gap].index:
                row[f"gap{gap}"] = None
                continue
            r = frames[gap].loc[read_id]
            row[f"gap{gap}"] = {"tl_bp": int(r.telomere_length), "read_bp": int(r.read_length),
                                   "mapq": int(r.mapping_quality), "chr": str(r.chromosome),
                                   "arm": str(r.arm), "direction": str(r.direction)}
        rows.append(row)
    comparisons = []
    for a, b in ((20, 100), (20, 250), (100, 250)):
        ia, ib = set(frames[a].index), set(frames[b].index)
        matched = sorted(ia & ib)
        delta = (frames[b].loc[matched, "telomere_length"] - frames[a].loc[matched, "telomere_length"]).to_numpy()
        coordinate_changed = sum(
            any(frames[a].loc[r, k] != frames[b].loc[r, k] for k in ("chromosome", "arm", "direction", "reference_start", "reference_end"))
            for r in matched)
        comparisons.append({"from_gap": a, "to_gap": b, "matched_n": len(matched),
                            "from_only_n": len(ia-ib), "to_only_n": len(ib-ia),
                            "changed_length_n": int((delta != 0).sum()),
                            "changed_alignment_n": coordinate_changed,
                            "median_delta_bp": float(np.median(delta)),
                            "max_abs_delta_bp": int(np.abs(delta).max()),
                            "mean_delta_bp": float(np.mean(delta))})
    result = {
        "schema_version": 1,
        "metadata": {
            "repository": REPO, "commit": COMMIT,
            "repository_license": "MIT (Copyright 2023 santiago-es)",
            "bam_sha256": BAM_SHA, "wheel_sha256": WHEEL_SHA, "wheel_url": WHEEL_URL,
            "versions": {name: importlib.metadata.version(name) for name in ("telometer", "numpy", "pandas", "pysam", "regex", "scipy")},
            "min_read_bp": 1000, "gap_values_bp": list(GAPS), "threads": 2,
            "enrichment": False,
            "note": "Current-version test of the public subsampled capture BAM; not 2024 paper reproduction, not an independent biological replicate, and not valid for telomeres-per-Gb enrichment.",
            "gap_note": "Version 2.0.3 CLI default is 20 bp. The pinned repository README separately mentions 100 bp in a workflow example and 250 bp in its options table. All three configurations were explicit, not silently inferred.",
            "comparison_note": "A changed result is a sensitivity flag, not evidence that one setting is correct. Telometer chooses one alignment per read ID; compare per-setting coordinate/arm fields because the selected alignment can also change."
        },
        "bam": bam_stats,
        "summaries": {str(g): summary(frames[g]) for g in GAPS},
        "comparisons": comparisons,
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("bam", "summaries", "comparisons")}, indent=2))


if __name__ == "__main__":
    main()
