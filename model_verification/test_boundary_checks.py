"""Independent finite-boundary and qualified-theory checks."""

from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import math
import random
import unittest

import numpy as np

import make_figures as model
import reproduce_numerical_checks as numerical


def decimal_rate(annual, interval):
    with localcontext() as context:
        context.prec = 110
        p, t = Decimal.from_float(annual), Decimal.from_float(interval)
        return float(1 - ((1 - p).ln() * t).exp())


class BoundaryCheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = {row["Id"]: row for row in numerical.native_checks(boundary=True)}

    def test_complete_native_regression_set(self):
        self.assertEqual(len(self.native), 21)
        self.assertTrue(all(row["Passed"] for row in self.native.values()))

    def test_exact_capacity_boundary_period(self):
        design = dict(mcov=1, macc=0, mtime=1, mplur=1, gamma=1, rho=1)
        parameters = dict(model.BASE, eta=.9, delta0=.1, lockR=0, chi=0, phi=0,
                          theta=1, q=1, sig=0, noise=1, floors=(0,) * 6,
                          salience_form="linear_capped")
        result = model.period(design, parameters, [.3] * 6, [.55, .6], 0)
        self.assertEqual(result["Cn"], [0.0] * 6)
        self.assertEqual(result["raw_capacity"], [0.0] * 6)
        self.assertEqual(result["Cn"], self.native["FR02_capacity_boundary_period"]["Values"])

    def test_exact_capacity_boundary_feedback(self):
        record = dict(
            closure="capacity_subsystem_fixed_forcing", eta=.9, delta0=.1, chi=0,
            lockR=0, theta=1, phi=0, mcov=1, macc=0, mtime=1, mplur=1, gamma=1,
            rho=1, h=2, s0=.4, kappa=.25, q=[1] * 6, beta=[1 / 6] * 6,
            external_capacity_change=[0] * 6, required_routes=[True], available_routes=[True],
            monitoring_form="product", salience_form="linear_capped", signal=[0], error=[1],
            AI=[[0]], AT=[[0]])
        np.testing.assert_array_equal(numerical.feedback_map(record, [.3] * 6), np.zeros(6))

    def test_nonnegative_update_against_decimal_oracle(self):
        rng = random.Random(20261024)
        with localcontext() as context:
            context.prec = 100
            for trial in range(1000):
                eta = rng.random()
                decay = (1 - eta) * rng.random()
                consumption = max(0, 1 - eta - decay) * rng.random()
                stock, retained, relevance, valid, accuracy = [rng.random() for _ in range(5)]
                result = model._capacity_coordinate(
                    stock, eta, retained, relevance, valid, accuracy, decay, consumption)
                c, a, r, q, v, accuracy_d, d, z = map(
                    Decimal.from_float, (stock, eta, retained, relevance, valid, accuracy, decay, consumption))
                expected = c + a * r * q * (v * (1 - c) - (1 - v) * (1 - accuracy_d) * c) - (d + z) * c
                with self.subTest(trial=trial):
                    self.assertGreaterEqual(result, 0)
                    self.assertLessEqual(result, 1)
                    self.assertLessEqual(abs(result - float(expected)), 4e-16)

    def test_tiny_annual_rate_long_interval(self):
        expected = decimal_rate(1e-18, 1e18)
        self.assertAlmostEqual(model.review_rate(1e-18, 1e18), expected, places=15)
        self.assertAlmostEqual(self.native["FR04_tiny_annual_rate"]["Values"][0], expected, places=15)

    def test_tiny_interval_nonzero_rate(self):
        expected = decimal_rate(.2, 1e-18)
        self.assertGreater(model.review_rate(.2, 1e-18), 0)
        self.assertLessEqual(abs(model.review_rate(.2, 1e-18) - expected), 1e-33)
        self.assertLessEqual(abs(self.native["FR05_tiny_review_interval"]["Values"][0] - expected), 1e-33)

    def test_rate_endpoints_and_one_year_identity(self):
        self.assertEqual(model.review_rate(0, 1e308), 0)
        self.assertEqual(model.review_rate(1, 1e-308), 1)
        for rate in (1e-300, .1, .9999999999999999):
            self.assertEqual(model.review_rate(rate, 1), rate)

    def test_rate_guards(self):
        for annual, interval in ((-.1, 1), (1.1, 1), (True, 1), (.1, 0),
                                 (.1, -1), (.1, math.inf), (math.nan, 1)):
            with self.subTest(annual=annual, interval=interval), self.assertRaises(ValueError):
                model.review_rate(annual, interval)

    def test_long_interval_actual_decay(self):
        parameters = dict(model.BASE, delta0=1e-18, period=1e18, eta=0, lockR=0,
                          chi=0, phi=0, sig=0, noise=0, floors=(0,) * 6)
        result = model.period(model.DESIGNS["B"], parameters, [.4] * 6, [.55, .6], 0)
        self.assertAlmostEqual(result["Cn"][0], .4 * (1 - decimal_rate(1e-18, 1e18)), places=15)

    def test_converted_rate_guard_cannot_be_bypassed(self):
        parameters = dict(model.BASE, delta0=1e-18, period=1e18, eta=.7, lockR=0,
                          chi=0, phi=0, sig=0, noise=0, floors=(0,) * 6)
        with self.assertRaisesRegex(ValueError, "S29"):
            model.period(model.DESIGNS["B"], parameters, [.4] * 6, [.55, .6], 0)

    def test_signed_fusion_permutations_and_duplicate_insertions(self):
        values, identities = [1e308, 1.0, -1e308], ["a:v1", "b:v1", "c:v1"]
        expected = float(sum(map(Fraction.from_float, values)) / 3)
        checked = 0
        for order in itertools.permutations(range(3)):
            candidates = [list(order)]
            for copied in range(3):
                for position in range(4):
                    records = list(order)
                    records.insert(position, copied)
                    candidates.append(records)
            for records in candidates:
                count = len(records)
                actual = model.canonical_evidence(
                    [values[i] for i in records], [1] * count, [0] * count,
                    [1] * count, [identities[i] for i in records])[0]
                self.assertEqual(actual, expected)
                checked += 1
        self.assertEqual(checked, 78)
        self.assertEqual(self.native["FR03_signed_fusion_invariance"]["Values"], [expected, 78])

    def test_certain_harm_rejects_every_subunit_tolerance(self):
        bounds = self.native["FR10_certain_harm"]["Values"][:3]
        self.assertEqual(bounds, [1, 1, 1])
        self.assertFalse(model.risk_accepts(bounds, [math.nextafter(1, 0)] * 3))

    def test_zero_probability_safe_scenario_preserves_certainty(self):
        self.assertEqual(self.native["FR11_zero_probability_safe_scenario"]["Values"], [1, 1, 1])
        self.assertEqual(self.native["FR12_zero_harm"]["Values"], [0, 0, 0])

    def test_hill_adjacent_floats_high_precision(self):
        scale, strength, exponent = 1e200, math.nextafter(1e200, math.inf), 1e16
        with localcontext() as context:
            context.prec = 110
            logarithm = (Decimal.from_float(strength) / Decimal.from_float(scale)).ln()
            z = Decimal.from_float(exponent) * logarithm
            expected = float(1 / (1 + (-z).exp()))
        self.assertLessEqual(abs(model.salience(strength, scale, exponent) - expected), 2e-15)
        self.assertLessEqual(abs(self.native["FR13_near_threshold_hill"]["Values"][0] - expected), 2e-15)

    def test_hill_exact_unit_rescaling(self):
        original = model.salience(math.nextafter(1024.0, math.inf), 1024.0, 1e16)
        self.assertEqual(original, model.salience(math.nextafter(1.0, math.inf), 1.0, 1e16))
        self.assertEqual(self.native["FR14_hill_unit_invariance"]["Values"][0],
                         self.native["FR14_hill_unit_invariance"]["Values"][1])

    def test_hill_wide_dynamic_range(self):
        self.assertEqual(model.salience(1e308, 1e-308, 1e308), 1)
        self.assertEqual(model.salience(1e-308, 1e308, 1e308), 0)
        self.assertEqual(model.salience(5e-324, 5e-324, 1e16), .5)

    def test_critical_floor_allows_above_floor_tradeoff_only(self):
        passes = lambda performance: model.protection_holds(
            performance, (.4, 0), (.5,) * 6, (0,) * 6, (), ())
        self.assertTrue(passes((.8, .4)))
        self.assertTrue(passes((.7, .5)))
        self.assertAlmostEqual(.5 * (.8 + .4), .5 * (.7 + .5))
        self.assertFalse(passes((.39, 1)))

    def test_error_counterexample_preserves_the_qualified_theorem(self):
        design = dict(mcov=1, macc=.8, mtime=1, mplur=1, gamma=1, rho=1)
        parameters = dict(model.BASE, eta=.5, delta0=0, lockR=0, chi=0, phi=0,
                          theta=1, q=1, sig=.1, noise=0, s0=1, floors=(0,) * 6)
        change = lambda p: model.period(design, p, [.1] * 6, [.55, .6], 1)["Cn"][0] - .1
        clean = change(parameters)
        noisy = change(dict(parameters, noise=.5))
        pure = change(dict(parameters, sig=0, noise=.5))
        self.assertAlmostEqual(clean, float(Fraction(9, 2020)), places=15)
        self.assertAlmostEqual(noisy, float(Fraction(13, 945)), places=15)
        self.assertAlmostEqual(pure, float(Fraction(-1, 500)), places=15)
        self.assertGreater(noisy, clean)

    def test_taxonomy_dependence_requires_allocation_assumptions(self):
        def threshold(beta, q):
            qv = math.fsum(b * relevance * .5 for b, relevance in zip(beta, q))
            return model.critical_learning(.5, qv, .01)["value"]
        six = threshold([1 / 6] * 6, [1, 0, 0, 0, 0, 0])
        seven = threshold([1 / 6, 1 / 12, 1 / 12, 1 / 6, 1 / 6, 1 / 6, 1 / 6],
                          [1, 0, 0, 0, 0, 0, 0])
        self.assertAlmostEqual(six, .24)
        self.assertAlmostEqual(seven, six)
        self.assertAlmostEqual(threshold([1 / 6] * 6, [1 / 6] * 6), .24)
        self.assertAlmostEqual(threshold([1 / 7] * 7, [1 / 7] * 7), .28)

    def test_zero_evidence_has_zero_validity_and_learning(self):
        parameters = dict(model.BASE, sig=0, noise=0, chi=0, lockR=0)
        result = model.period(model.DESIGNS["B"], parameters, [.4] * 6, [.55, .6], 0)
        self.assertEqual(result["v"], 0)
        self.assertEqual(result["ell"], 0)
        retained = model.DESIGNS["B"]["rho"] * (parameters["theta"] + (1 - parameters["theta"]) * .4)
        decay = parameters["delta0"] * (1 - parameters["phi"] * retained)
        self.assertAlmostEqual(result["Cn"][0], .4 * (1 - decay), places=15)


if __name__ == "__main__":
    unittest.main()
