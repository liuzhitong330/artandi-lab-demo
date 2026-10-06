#!/usr/bin/env python3
"""Independent standard-library checks of the derived Telometer audit JSON."""
import json
import math
from pathlib import Path
import statistics

path = Path(__file__).resolve().parents[1] / "data" / "telometer-audit.json"
d = json.loads(path.read_text())
assert len(d["rows"]) == 738
assert len({r["read_id"] for r in d["rows"]}) == 738
assert len({r["id"] for r in d["rows"]}) == 738
assert d["bam"]["alignment_records"] == 1640
assert d["bam"]["distinct_read_ids"] == 812
assert d["metadata"]["enrichment"] is False
assert d["metadata"]["versions"]["telometer"] == "2.0.3"

for gap, count, median in ((20, 738, 3876), (100, 736, 3900), (250, 733, 3924)):
    values = [r[f"gap{gap}"]["tl_bp"] for r in d["rows"] if r[f"gap{gap}"] is not None]
    summary = d["summaries"][str(gap)]
    assert len(values) == count == summary["n"]
    assert statistics.median(values) == median == summary["median_bp"]
    assert math.isclose(statistics.mean(values), summary["mean_bp"], rel_tol=1e-12)
    assert sum(x < 2000 for x in values) == summary["below_2000_n"]
    assert min(values) == summary["min_bp"] and max(values) == summary["max_bp"]
    for row in d["rows"]:
        value = row[f"gap{gap}"]
        if value is not None:
            assert 0 < value["tl_bp"] < value["read_bp"]
            assert 0 <= value["mapq"] <= 255

for comparison in d["comparisons"]:
    a, b = comparison["from_gap"], comparison["to_gap"]
    av = {r["read_id"]: r[f"gap{a}"] for r in d["rows"] if r[f"gap{a}"] is not None}
    bv = {r["read_id"]: r[f"gap{b}"] for r in d["rows"] if r[f"gap{b}"] is not None}
    common = set(av) & set(bv)
    delta = [bv[r]["tl_bp"] - av[r]["tl_bp"] for r in common]
    assert len(common) == comparison["matched_n"]
    assert len(set(av) - set(bv)) == comparison["from_only_n"]
    assert len(set(bv) - set(av)) == comparison["to_only_n"]
    assert sum(x != 0 for x in delta) == comparison["changed_length_n"]
    assert statistics.median(delta) == comparison["median_delta_bp"]
    assert max(map(abs, delta)) == comparison["max_abs_delta_bp"]
    assert math.isclose(statistics.mean(delta), comparison["mean_delta_bp"], rel_tol=1e-12)

delta_100_250 = [r["gap250"]["tl_bp"] - r["gap100"]["tl_bp"] for r in d["rows"]
                 if r["gap100"] is not None and r["gap250"] is not None]
assert sum(abs(x) >= 100 for x in delta_100_250) == 243
assert sum(abs(x) >= 500 for x in delta_100_250) == 11
assert sum(abs(x) >= 1000 for x in delta_100_250) == 3
print("Telometer audit: independent summary, null-handling and sensitivity checks passed.")
