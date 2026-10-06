"""Offline regression and boundary checks for the published row-level dataset."""
import csv
import importlib.util
from itertools import combinations
import json
import math
from pathlib import Path
import statistics
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_tertko", ROOT / "scripts" / "build_tertko.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)
DATA = json.loads((ROOT / "data" / "tertko.json").read_text())
CHECKS = json.loads((ROOT / "data" / "analysis_checks.json").read_text())


class TertkoTests(unittest.TestCase):
    def test_row_identity_and_counts(self):
        self.assertEqual(DATA["columns"], ["row", "day", "telomere_bp", "read_bp", "chromosome", "arm", "direction"])
        self.assertEqual(len(DATA["rows"]), 2610)
        self.assertEqual([r[0] for r in DATA["rows"]], list(range(3, 2613)))
        self.assertEqual(DATA["rows"][0], [3, 66, 3415, 7612, "chr1", "p", "fwd"])
        self.assertEqual(DATA["rows"][-1], [2612, 98, 1228, 9981, "altX_q_0", "q", "fwd"])
        self.assertEqual(DATA["days"], [66, 78, 98, 105])
        self.assertEqual({day: sum(r[1] == day for r in DATA["rows"]) for day in DATA["days"]},
                         {66: 377, 78: 350, 98: 968, 105: 915})

    def test_numeric_integrity_and_flags(self):
        for row in DATA["rows"]:
            self.assertEqual(len(row), 7)
            self.assertTrue(all(v is not None for v in row))
            self.assertTrue(all(math.isfinite(row[i]) and row[i] > 0 for i in [2, 3]))
            self.assertLessEqual(row[2], row[3])
        self.assertEqual(CHECKS["audit"]["nonpositive_or_nonfinite_length_rows"], [])
        self.assertEqual(CHECKS["audit"]["telomere_exceeds_read_length_rows"], [])

    def test_repeat_rows_preserved_not_resolved(self):
        visible = [tuple(r[1:]) for r in DATA["rows"]]
        self.assertEqual(len(visible) - len(set(visible)), 25)
        self.assertEqual(DATA["duplicate_visible_rows"], 25)
        self.assertEqual(len(CHECKS["audit"]["repeated_visible_row_groups"]), 25)

    def test_summary_against_independent_aggregates(self):
        # Independently aggregate JSON rows using the standard library rather
        # than comparing the build function's output to itself.
        known = {
            66: (377, 4296.909814323608, 3886, 32),
            78: (350, 3345.6085714285714, 2824.5, 95),
            98: (968, 2617.240702479339, 2251, 409),
            105: (915, 2639.3857923497267, 2066, 423),
        }
        for summary in DATA["summary"]:
            day = summary["day"]
            values = [r[2] for r in DATA["rows"] if r[1] == day]
            n, mean, median, short_n = known[day]
            self.assertEqual(summary["n"], n)
            self.assertAlmostEqual(summary["mean"], sum(values) / len(values), places=10)
            self.assertAlmostEqual(summary["mean"], mean, places=10)
            self.assertEqual(summary["median"], statistics.median(values))
            self.assertEqual(summary["median"], median)
            self.assertEqual(summary["short_lt_2000_n"], short_n)
            self.assertEqual(summary["short_lt_2000_fraction"], short_n / n)

    def test_quantile_and_strict_boundary_contract(self):
        self.assertIsNone(builder.quantile([], 0.5))
        self.assertEqual(builder.quantile([5], 0.25), 5)
        self.assertEqual(builder.quantile([0, 10, 20, 30], 0.25), 7.5)
        synthetic = [[3, 66, 1999, 4000, "chr1", "p", "fwd"],
                     [4, 66, 2000, 3999, "chr1", "p", "fwd"],
                     [5, 66, 2001, 4000, "chr1", "p", "fwd"]]
        all_rows = builder.describe(synthetic, 66, threshold=2000)
        self.assertEqual(all_rows["short_lt_2000_n"], 1)
        retained = builder.describe(synthetic, 66, threshold=2000, min_read=4000)
        self.assertEqual(retained["n"], 2)
        self.assertEqual(retained["short_lt_2000_fraction"], 0.5)
        empty = builder.describe(synthetic, 105)
        self.assertEqual(empty["n"], 0)
        self.assertIsNone(empty["mean"])
        self.assertIsNone(empty["short_lt_2000_fraction"])

    def test_all_sensitivity_cases(self):
        expected = {(a, b, threshold, minimum)
                    for a, b in combinations(DATA["days"], 2)
                    for threshold in [1000, 1500, 2000, 2500, 3000, 4000]
                    for minimum in [0, 1000, 4000, 8000]}
        observed = set()
        for case in CHECKS["contrasts"]:
            a, b, threshold, minimum = (case[k] for k in
                ["earlier_day", "later_day", "threshold_bp", "minimum_read_bp"])
            observed.add((a, b, threshold, minimum))
            earlier = [r[2] for r in DATA["rows"] if r[1] == a and r[3] >= minimum]
            later = [r[2] for r in DATA["rows"] if r[1] == b and r[3] >= minimum]
            self.assertEqual(case["earlier_n"], len(earlier))
            self.assertEqual(case["later_n"], len(later))
            self.assertAlmostEqual(case["mean_change_bp"], sum(later) / len(later) - sum(earlier) / len(earlier), places=10)
            self.assertAlmostEqual(case["median_change_bp"], statistics.median(later) - statistics.median(earlier), places=10)
            delta = 100 * (sum(v < threshold for v in later) / len(later) - sum(v < threshold for v in earlier) / len(earlier))
            self.assertAlmostEqual(case["short_fraction_change_pp"], delta, places=10)
        self.assertEqual(observed, expected)
        self.assertEqual(len(CHECKS["contrasts"]), len(expected))

    def test_highlighted_endpoint_disagreement(self):
        case = next(c for c in CHECKS["contrasts"] if
                    (c["earlier_day"], c["later_day"], c["threshold_bp"], c["minimum_read_bp"]) == (98, 105, 2000, 0))
        self.assertAlmostEqual(case["mean_change_bp"], 22.14508987038787)
        self.assertEqual(case["median_change_bp"], -185)
        self.assertAlmostEqual(case["short_fraction_change_pp"], 3.977442081018835)
        self.assertTrue(case["mean_median_opposite_signs"])

    def test_export_and_provenance(self):
        with (ROOT / "data" / "tertko.csv").open(newline="") as stream:
            csv_rows = list(csv.reader(stream))
        self.assertEqual(csv_rows[0], DATA["columns"])
        self.assertEqual(csv_rows[1:], [[str(v) for v in row] for row in DATA["rows"]])
        self.assertEqual(DATA["source"]["sha256"], builder.WORKBOOK_SHA256)
        self.assertEqual(DATA["source"]["archive_sha256"], builder.SOURCE_ARCHIVE_SHA256)
        self.assertEqual(DATA["source"]["sheet"], "hesc_tertko")
        self.assertIn("105", CHECKS["audit"]["source_day_discrepancy"])
        self.assertIn("108", CHECKS["audit"]["source_day_discrepancy"])


if __name__ == "__main__":
    unittest.main()
