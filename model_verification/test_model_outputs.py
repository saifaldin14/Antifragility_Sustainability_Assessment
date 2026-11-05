import csv
import json
import math
from pathlib import Path
import tempfile
import unittest

import make_figures as model
import reproduce_scenario_analyses as scenarios


class ModelTests(unittest.TestCase):
    def test_preserved_annual_values(self):
        model.check_against_script()

    def test_preserved_latent_classes(self):
        for design, expected in (("A", "UUXXX"), ("B", "GGGGG")):
            result = model.run_design(design, dict(model.BASE), 5, lambda t: True)
            self.assertEqual(result["classes"], expected)
            self.assertEqual(result["evidence_status"], "NotAssessed_IllustrativeSimulation")

    def test_cap_comparison(self):
        expected = {
            (1.5, "A"): (0.4773, 0.3627),
            (1.5, "B"): (0.5345, 0.4013),
            (math.inf, "A"): (0.3525, 0.3627),
            (math.inf, "B"): (0.4099, 0.4066),
        }
        for (cap, design), (social, capacity) in expected.items():
            with self.subTest(cap=cap, design=design):
                result = model.period(model.DESIGNS[design], dict(model.BASE, period=5, cap=cap),
                                      [0.4] * 6, [0.55, 0.60], 1)
                self.assertAlmostEqual(result["Sn"][0], social, places=4)
                self.assertAlmostEqual(sum(result["Cn"]) / 6, capacity, places=4)
                self.assertFalse(result["SB"])
                self.assertEqual(result["d"], 1.5 if cap == 1.5 else 4)

    def test_zero_disturbance_counterfactual(self):
        result = model.period(model.DESIGNS["B"], dict(model.BASE, period=5, cap=math.inf),
                              [0.4] * 6, [0.55, 0.60], 0)
        self.assertEqual(result["d"], 0)

    def test_invalid_caps(self):
        for cap in (0, -1, -math.inf, math.nan):
            with self.subTest(cap=cap), self.assertRaises(ValueError):
                model.period(model.DESIGNS["B"], dict(model.BASE, cap=cap), [0.4] * 6, [0.55, 0.60], 1)

    def test_baseline_is_outside_low_strength_domain(self):
        self.assertGreater(math.hypot(0.8, 0.2), 0.4 / math.sqrt(3))

    def test_export_shapes_and_status(self):
        previous = model.OUT
        try:
            with tempfile.TemporaryDirectory() as directory:
                model.OUT = Path(directory)
                model.export_model_records()
                for name, count in (("Cap_Comparison.csv", 4), ("Illustrative_Trajectories.csv", 62)):
                    with (model.OUT / name).open(newline="") as handle:
                        rows = list(csv.DictReader(handle))
                    self.assertEqual(len(rows), count)
                    self.assertTrue(all(row["evidence_status"] == model.SIMULATION_EVIDENCE_STATUS for row in rows))
        finally:
            model.OUT = previous

    def test_import_does_not_use_test_runner_arguments(self):
        self.assertEqual(model.OUT, Path("."))

    def test_routine_learning_without_disturbance_has_no_attribution(self):
        c, s = [0.4] * 6, [0.55, 0.60]
        p = dict(model.BASE, routine=(0.2, 0.0), noise=0.0)
        actual = model.period(model.DESIGNS["B"], p, c, s, 0)
        counter = model.period(model.DESIGNS["B"], p, c, s, 0)
        counter["C"] = c
        self.assertGreater(sum(actual["Cn"]), sum(c))
        self.assertEqual(actual["Cn"], counter["Cn"])
        self.assertEqual(model.latent_class_of(actual, counter), "R")

    def test_routine_aware_contrasts(self):
        rows = scenarios.structural_rows()[:4]
        self.assertTrue(all(a["contrast5"] > b["contrast5"] for a, b in zip(rows, rows[1:])))
        self.assertAlmostEqual(rows[1]["capacity5"], 0.4958014119755938, places=12)
        self.assertAlmostEqual(rows[1]["comparator5"], 0.39255867919684384, places=12)

    def test_component_histories_are_complete(self):
        result = model.run_design("B", dict(model.BASE), 5, lambda _: True)
        for key in ("components", "components_cf"):
            self.assertEqual(len(result[key]), 6)
            self.assertTrue(all(len(row) == 6 for row in result[key]))
        self.assertAlmostEqual(result["components"][-1][model.iR], 0.44877296084804047, places=12)

    def test_component_floor_rejects_without_clipping(self):
        p = dict(model.BASE, floors=(0, 0, 0, 0, 0.6, 0))
        result = model.run_design("B", p, 5, lambda _: True)
        self.assertEqual(result["classes"], "XXXXX")
        self.assertFalse(result["retained_gain"])
        self.assertLess(result["C"][model.iR], 0.6)

    def test_above_floor_tradeoff_is_explicitly_permitted(self):
        final = [0.6, 0.6, 0.6, 0.6, 0.3, 0.6]
        self.assertGreater(sum(final) / 6, 0.5)
        self.assertTrue(model.protection_holds((0.6,), (0.4,), final, (0, 0, 0, 0, 0.25, 0), (), ()))
        self.assertFalse(model.protection_holds((0.6,), (0.4,), final, (0, 0, 0, 0, 0.4, 0), (), ()))

    def test_episodic_horizon_fails_protection(self):
        result = model.run_design("B", dict(model.BASE), 30, lambda t: t % 5 == 2)
        self.assertEqual(result["classes"], "RRGRRRRGRRRRGRRRXGRXXXXXXXXXXX")
        self.assertEqual(result["classes"].count("G"), 4)
        self.assertFalse(result["retained_gain"])

    def test_mature_episodic_start_is_not_protected(self):
        p, event = dict(model.BASE), lambda t: t % 5 == 2
        initial = model.run_design("B", p, 4000, event)["C"]
        self.assertLess(initial[model.iR], 0.3)
        result = model.run_design("B", p, 30, event, C0=initial)
        self.assertEqual(result["classes"], "X" * 30)
        self.assertFalse(result["protected_horizon"])
        self.assertFalse(result["retained_gain"])

    def test_canonical_units_and_redundancy(self):
        self.assertEqual(model.canonical_evidence([20], [100], [0], [1], ["a:v1"]), (0.2,))
        self.assertEqual(model.canonical_evidence([2000], [10000], [0], [1], ["a:v1"]), (0.2,))
        self.assertEqual(model.canonical_evidence([20, 20], [100, 100], [0, 0], [1, 1], ["a:v1", "a:v1"]), (0.2,))

    def test_canonical_input_guards(self):
        cases = [
            ([], [], [], []), ([1], [0], [0], [1]), ([math.nan], [1], [0], [1]),
            ([1], [1], [1], [1]), ([1], [1], [10**9], [1]),
            ([1], [1], [0], [0]), ([1], [1], [0], [-1]),
            ([1], [1], [True], [1]), ([1], [1], [0], [math.inf]),
            ([1, 2], [1], [0], [1]),
        ]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                model.canonical_evidence(*args, ["a:v1"] * len(args[0]))

    def test_monitoring_forms_keep_missing_channel_veto(self):
        for form in ("product", "geometric", "minimum"):
            self.assertEqual(model.monitoring_quality((0.9, 0, 0.7), form), 0)
        with self.assertRaises(ValueError):
            model.monitoring_quality((0.9, 0.8, 0.7), "unknown")

    def test_temporary_recovery_does_not_erase_breach(self):
        self.assertFalse(model.protection_holds((0.6,), (0.4,), (0.5,), (0.3,), (), (),
                                               peak=(0.3,), peak_floors=(0.4,)))

    def test_cumulative_budget_is_binding(self):
        self.assertFalse(model.protection_holds((0.6,), (0.4,), (0.5,), (0.3,), (), (),
                                               cumulative=(1.2,), cumulative_limits=(1.0,)))

    def test_invalid_protection_vectors(self):
        with self.assertRaises(ValueError):
            model.protection_holds((0.6,), (0.4,), (0.5,), (), (), ())
        with self.assertRaises(ValueError):
            model.protection_holds((0.6,), (0.4,), (math.nan,), (0.3,), (), ())

    def test_structural_resource_boundary(self):
        rows = {row["scenario"]: row for row in scenarios.structural_rows()}
        self.assertEqual(len(rows), 13)
        for name in ("coupled_resource_half", "coupled_resource_quarter"):
            self.assertEqual(rows[name]["classes"], "RRRRR")
            self.assertGreater(rows[name]["contrast5"], 0)
            self.assertLess(rows[name]["capacity5"], 0.4)

    def test_matched_benchmarks(self):
        rows = scenarios.matched_benchmarks()
        self.assertEqual(len(rows), 5)
        self.assertEqual([row["verdict"] for row in rows],
                         ["Robust", "CapabilityGain", "Robust", "Unacceptable", "Unacceptable"])
        self.assertTrue(all(row["resource_units"] == 1 for row in rows))

    def test_completed_dossier_rejects_positive_gain(self):
        rows = scenarios.dossier_rows()
        self.assertEqual(len(rows), 5)
        self.assertGreater(rows[-1]["capacity_index"], rows[-1]["comparator_index"])
        self.assertAlmostEqual(rows[-1]["cumulative_burden"], 1.2)
        self.assertFalse(rows[-1]["all_declared_protections"])
        self.assertEqual(rows[3]["instruction"], "Escalate before next planned increment")

    def test_native_scenario_checks_and_exact_vectors(self):
        self.assertEqual(len(scenarios.native_checks()), 31)

    def test_all_native_draws_cross_implementation(self):
        rows = scenarios.extract_draws(scenarios.read_log(
            Path(__file__).parent / "Verified_Native_Output.txt"))
        checked = scenarios.check_global(rows)
        self.assertEqual(checked["numeric_comparisons"], 8000)
        self.assertEqual(checked["class_sequence_comparisons"], 4000)
        self.assertEqual(checked["counts"]["BcapabilityGainEveryYear"], 358)
        self.assertEqual(checked["counts"]["BcapabilityGainAnyYear"], 447)
        self.assertEqual(checked["percentages"]["BaboveCounterfactual"], 97.1)
        self.assertLessEqual(checked["maximum_absolute_difference"], 1e-12)

    def test_tampered_draw_is_detected(self):
        rows = scenarios.extract_draws(scenarios.read_log(
            Path(__file__).parent / "Verified_Native_Output.txt"))
        rows[0]["B_capacity5"] = str(float(rows[0]["B_capacity5"]) + 0.01)
        with self.assertRaises(AssertionError):
            scenarios.check_global(rows)

    def test_incomplete_draws_are_rejected(self):
        with self.assertRaises(ValueError):
            scenarios.extract_draws("GD|draw\nGD|1\n")
        with self.assertRaises(ValueError):
            scenarios.check_global([])

    def test_utf16_log_is_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "native.txt"
            path.write_text("RESULT: PASS", encoding="utf-16")
            self.assertEqual(scenarios.read_log(path), "RESULT: PASS")

    def test_generated_scenario_record_shapes(self):
        root = Path(__file__).parent
        for filename, expected in (("Global_Sensitivity_Draws.csv", 2000),
                                   ("Structural_Robustness.csv", 13),
                                   ("Matched_Decision_Benchmarks.csv", 5),
                                   ("Illustrative_Dossier.csv", 5)):
            with (root / filename).open(encoding="utf-8", newline="") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), expected)
        cases = json.loads((root / "Scenario_Checks.json").read_text(encoding="utf-8"))
        self.assertEqual(len(cases), 31)
        self.assertTrue(all(row["Passed"] for row in cases))


if __name__ == "__main__":
    unittest.main()
