import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np

import reproduce_cases as model


ROOT = Path(__file__).resolve().parent


class StatisticsTests(unittest.TestCase):
    def test_perfect_agreement(self):
        for metric in (model.ordinal_alpha, model.weighted_kappa):
            self.assertEqual(metric([1, 2, 3, 6], [1, 2, 3, 6]), 1)

    def test_constant_agreement_is_undefined(self):
        self.assertIsNone(model.ordinal_alpha([4, 4], [4, 4]))
        self.assertIsNone(model.weighted_kappa([4, 4], [4, 4]))

    def test_alpha_symmetry(self):
        a, b = [1, 3, 4, 6], [2, 3, 5, 4]
        self.assertEqual(model.ordinal_alpha(a, b), model.ordinal_alpha(b, a))

    def test_invalid_ratings(self):
        for values in ([], [0], [7], [2.5], [np.nan], [np.inf], [[1, 2]]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                model.ratings(values)

    def test_mismatched_pairs(self):
        with self.assertRaises(ValueError):
            model.contingency([1], [1, 2])

    def test_invalid_weight(self):
        with self.assertRaises(ValueError):
            model.weighted_kappa([1, 2], [1, 2], 3)

    def test_exact_requires_blocks(self):
        with self.assertRaises(ValueError):
            model.bootstrap([1, 2], [1, 2], model.ordinal_alpha, exact=True)

    def test_invalid_bootstrap_partition(self):
        with self.assertRaises(ValueError):
            model.bootstrap([1, 2], [1, 2], model.ordinal_alpha, [np.array([0])])


class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = model.read_csv(ROOT / "data" / "Ratings.csv")
        cls.evidence = model.read_csv(ROOT / "data" / "Case_Criteria.csv")
        cls.sources = model.read_csv(ROOT / "data" / "Source_Ledger.csv")
        with tempfile.TemporaryDirectory() as directory:
            cls.result = model.reproduce(ROOT / "data", Path(directory))

    def test_final_coefficients(self):
        result = self.result["final_round"]
        self.assertAlmostEqual(result["ordinal_alpha"], 0.443863770944856)
        self.assertAlmostEqual(result["quadratic_kappa"], 0.511960406928788)
        self.assertAlmostEqual(result["linear_kappa"], 0.314214463840399)

    def test_agreements_and_bias(self):
        result = self.result["final_round"]
        self.assertEqual(result["exact_agreements"], 28)
        self.assertEqual(result["within_one_category"], 67)
        self.assertAlmostEqual(result["mean_second_minus_first"], -0.386666666666667)

    def test_exact_case_interval(self):
        result = self.result["final_round"]["exact_case_bootstrap_alpha"]
        self.assertEqual(result["defined_resamples"] + result["undefined_resamples"], 3125)
        self.assertAlmostEqual(result["percentile_95"][0], -0.2160702980)
        self.assertAlmostEqual(result["percentile_95"][1], 0.6795442723)

    def test_revision_provenance(self):
        self.assertEqual(self.result["second_rater_revisions"], 30)
        self.assertEqual(self.result["revisions_upward"], 29)
        self.assertEqual(self.result["revisions_toward_first"], 21)
        self.assertEqual(self.result["recorded_consensus_cells"], 0)

    def test_baseline_reductions(self):
        comparison = model.compare_baseline(self.rows, model.read_csv(ROOT / "data" / "Thesis_Table20_Ratings.csv"))
        self.assertEqual(sum(row["change"] == -1 for row in comparison), 6)
        self.assertEqual(sum(row["change"] == 0 for row in comparison), 69)

    def test_no_case_causal_validation(self):
        result = model.assess_evidence(self.evidence, self.sources)
        self.assertTrue(all(row["inference_status"] == "InsufficientForCausalTest" for row in result.values()))

    def test_all_established_is_not_verified_effect(self):
        rows = copy.deepcopy(self.evidence)
        for row in rows:
            row["status"] = "E"
        result = model.assess_evidence(rows, self.sources)
        self.assertTrue(all(row["inference_status"] == "ReadyForIndependentEffectAndSafeguardAnalysis" for row in result.values()))

    def test_unknown_status_rejected(self):
        rows = copy.deepcopy(self.evidence)
        rows[0]["status"] = "?"
        with self.assertRaises(ValueError):
            model.assess_evidence(rows, self.sources)

    def test_duplicate_evidence_rejected(self):
        with self.assertRaises(ValueError):
            model.assess_evidence(self.evidence + [self.evidence[0]], self.sources)

    def test_missing_cell_rejected(self):
        with self.assertRaises(ValueError):
            model.assess_evidence(self.evidence[:-1], self.sources)

    def test_unknown_source_rejected(self):
        rows = copy.deepcopy(self.evidence)
        rows[0]["source_ids"] = "UNKNOWN"
        with self.assertRaises(ValueError):
            model.assess_evidence(rows, self.sources)

    def test_cross_case_source_rejected(self):
        rows = copy.deepcopy(self.evidence)
        rows[0]["source_ids"] = "A1"
        with self.assertRaises(ValueError):
            model.assess_evidence(rows, self.sources)

    def test_missing_limitation_rejected(self):
        rows = copy.deepcopy(self.evidence)
        rows[0]["limitation"] = ""
        with self.assertRaises(ValueError):
            model.assess_evidence(rows, self.sources)

    def test_duplicate_rating_rejected(self):
        with self.assertRaises(ValueError):
            model.load_ratings(self.rows + [self.rows[0]])

    def test_invented_consensus_rejected(self):
        rows = copy.deepcopy(self.rows)
        rows[0]["consensus"] = "4"
        with self.assertRaises(ValueError):
            model.load_ratings(rows)

    def test_ahrens_second_rater_has_one_level(self):
        self.assertEqual(self.result["final_round"]["by_case"]["C2"]["second_rater_distinct_levels"], 1)


if __name__ == "__main__":
    unittest.main()
