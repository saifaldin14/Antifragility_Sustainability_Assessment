import hashlib
import math
import unittest

import numpy as np

import make_figures as model
import reproduce_numerical_checks as numerical
import reproduce_scenario_analyses as verification


class NumericalCheckTests(unittest.TestCase):
    def test_native_source_provenance_accepts_exact_source(self):
        current = hashlib.sha256(verification.NATIVE.read_bytes()).hexdigest()
        verification.require_native_provenance("NATIVE_SOURCE_SHA256|" + current)

    def test_native_source_provenance_rejects_missing_stale_duplicate(self):
        current = hashlib.sha256(verification.NATIVE.read_bytes()).hexdigest()
        for text in ("", "NATIVE_SOURCE_SHA256|" + "0" * 64,
                     ("NATIVE_SOURCE_SHA256|" + current + "\n") * 2):
            with self.subTest(log=text):
                with self.assertRaises(ValueError):
                    verification.require_native_provenance(text)

    def test_exported_feedback_contract_and_boolean_state_rejection(self):
        record = dict(closure="capacity_subsystem_fixed_forcing", eta=.8, delta0=.08, chi=0,
                      lockR=0, theta=0, phi=0, mcov=1, macc=1, mtime=1, mplur=1, gamma=1,
                      rho=1, h=2, s0=.5, kappa=.2, q=[1] * 6, beta=[1/6] * 6,
                      external_capacity_change=[0] * 6, required_routes=[True], available_routes=[True],
                      monitoring_form="product", salience_form="hill", signal=[1], error=[0],
                      AI=[[0]], AT=[[0]])
        self.assertAlmostEqual(numerical.feedback_map(record, [.4] * 6)[0], .42944, places=14)
        self.assertAlmostEqual(numerical.feedback_map(dict(record, available_routes=[False]), [.4] * 6)[0],
                               .368, places=14)
        with self.assertRaises(ValueError):
            numerical.feedback_map(record, [True] * 6)

    def test_mixed_source_copy_invariance(self):
        original = model.canonical_evidence([.2, .6], [1, 1], [0, 0], [1, 1], ["a:v1", "b:v1"])
        copied = model.canonical_evidence([.2, .6, .2], [1, 1, 1], [0, 0, 0], [1, 1, 1],
                                          ["a:v1", "b:v1", "a:v1"])
        self.assertEqual(original, copied)
        self.assertEqual(original, (.4,))
        self.assertEqual(model.salience(model.stable_norm(original), .4, 2), .5)

    def test_frozen_source_weights(self):
        expected = model.canonical_evidence([.2, .6], [1, 1], [0, 0], [3, 1], ["a:v1", "b:v1"])
        actual = model.canonical_evidence([.2, .6, .2, .2], [1] * 4, [0] * 4, [3, 1, 3, 3],
                                          ["a:v1", "b:v1", "a:v1", "a:v1"])
        self.assertEqual(expected, actual)
        self.assertAlmostEqual(actual[0], .3)

    def test_source_conflicting_values_or_weights_rejected(self):
        for values, weights in (([.2, .6], [1, 1]), ([.2, .2], [1, 2])):
            with self.subTest(values=values, weights=weights), self.assertRaises(ValueError):
                model.canonical_evidence(values, [1, 1], [0, 0], weights, ["a:v1", "a:v1"])

    def test_distinct_sources_are_not_discarded(self):
        result = model.canonical_evidence([.2, .2, .6], [1] * 3, [0] * 3, [1] * 3,
                                          ["a:v1", "c:v1", "b:v1"])
        self.assertAlmostEqual(result[0], 1 / 3)

    def test_mixed_channels_and_unit_rescaling(self):
        original = model.canonical_evidence([20, 60, 30], [100, 100, 100], [0, 0, 1], [1] * 3,
                                            ["a:v1", "b:v1", "a:v1"])
        copied = model.canonical_evidence([2000, 6000, 30, 2000], [10000, 10000, 100, 10000],
                                          [0, 0, 1, 0], [1] * 4, ["a:v1", "b:v1", "a:v1", "a:v1"])
        self.assertEqual(original, copied)
        self.assertEqual(original, (.4, .3))

    def test_mixed_source_copy_downstream_invariance(self):
        original = model.canonical_evidence([.2, .6], [1, 1], [0, 0], [1, 1], ["a:v1", "b:v1"])
        copied = model.canonical_evidence([.2, .6, .2], [1] * 3, [0] * 3, [1] * 3,
                                          ["a:v1", "b:v1", "a:v1"])
        first = model.run_design("B", dict(model.BASE, routine=(original[0], 0)), 5, lambda _: True)
        second = model.run_design("B", dict(model.BASE, routine=(copied[0], 0)), 5, lambda _: True)
        for field in ("components", "components_cf", "performance", "classes", "retained_gain"):
            self.assertEqual(first[field], second[field])

    def test_large_source_weights(self):
        self.assertEqual(model.canonical_evidence([.2, .6], [1, 1], [0, 0], [1e308, 1e308],
                                                  ["a:v1", "b:v1"]), (.4,))

    def test_missing_source_identity_rejected(self):
        for identities in (None, [], [""], [" "]):
            with self.subTest(identities=identities), self.assertRaises(ValueError):
                model.canonical_evidence([.2], [1], [0], [1], identities)

    def test_stable_large_and_tiny_norm(self):
        self.assertEqual(model.stable_norm([1e200]), 1e200)
        self.assertEqual(model.stable_norm([1e-200]), 1e-200)
        self.assertEqual(model.stable_norm([]), 0)

    def test_norm_unrepresentable_result_is_explicit(self):
        with self.assertRaises(ValueError):
            model.stable_norm([1.7e308, 1.7e308])

    def test_stable_hill_extremes(self):
        self.assertEqual(model.salience(1e200, 1, 2), 1)
        self.assertEqual(model.salience(1e-200, 1e-200, 2), .5)
        self.assertEqual(model.salience(0, .4, 2), 0)
        self.assertEqual(model.salience(1e308, 1e-308, 1e308), 1)

    def test_stable_validity_extremes(self):
        self.assertEqual(model.validity([1e308], [1e308]), .5)
        self.assertEqual(model.validity([1e-308], [1e-308]), .5)
        self.assertEqual(model.validity([0], [0]), 0)

    def test_kernel_domain_errors(self):
        calls = (lambda: model.salience(-1, .4, 2), lambda: model.salience(1, 0, 2),
                 lambda: model.salience(1, .4, .5), lambda: model.validity([1], [1, 2]),
                 lambda: model.clip(math.nan), lambda: model.stable_norm([math.inf]))
        for call in calls:
            with self.subTest(call=call), self.assertRaises(ValueError):
                call()

    def test_negative_eta_is_rejected(self):
        with self.assertRaises(ValueError):
            model.period(model.DESIGNS["B"], dict(model.BASE, eta=-.1), [.4] * 6, [.55, .6], 1)

    def test_period_parameter_guards(self):
        invalid = dict(eta=1.1, delta0=-.01, theta=1.1, phi=-.1, chi=math.nan, lockR=1.1,
                       q=-.1, c0=1.1, h=.5, s0=0, kappa=0, period=0, sig=-1,
                       noise=math.inf, routine=(math.nan, 0), floors=(0,) * 5,
                       horizon_tolerance=-.1, required_routes=None, available_routes=(1,))
        for name, value in invalid.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                model.period(model.DESIGNS["B"], dict(model.BASE, **{name: value}), [.4] * 6, [.55, .6], 1)

    def test_governance_guards(self):
        for name in model.DESIGNS["B"]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                model.period(dict(model.DESIGNS["B"], **{name: math.nan}), model.BASE, [.4] * 6, [.55, .6], 1)

    def test_state_and_shape_guards(self):
        for c, s in (([.4] * 5, [.55, .6]), ([math.nan] * 6, [.55, .6]),
                     ([1.1] * 6, [.55, .6]), ([.4] * 6, [.55]), ([.4] * 6, [.55, -.1])):
            with self.subTest(c=c, s=s), self.assertRaises(ValueError):
                model.period(model.DESIGNS["B"], model.BASE, c, s, 1)

    def test_post_conversion_s29(self):
        with self.assertRaises(ValueError):
            model.period(model.DESIGNS["B"], dict(model.BASE, eta=.95, period=5), [.4] * 6, [.55, .6], 1)

    def test_classifier_rejects_nonfinite_before_label(self):
        step = model.period(model.DESIGNS["B"], model.BASE, [.4] * 6, [.55, .6], 1)
        counter = model.period(model.DESIGNS["B"], model.BASE, [.4] * 6, [.55, .6], 0)
        counter["C"] = [.4] * 6
        for name in ("MQ", "SC", "SCn", "SCmin", "d"):
            for invalid in (math.nan, math.inf, -math.inf):
                with self.subTest(name=name, invalid=invalid), self.assertRaises(ValueError):
                    model.latent_class_of(dict(step, **{name: invalid}), counter)

    def test_risk_tolerance_one_rejected(self):
        with self.assertRaises(ValueError):
            model.risk_accepts([1, 1, 1], [1, 1, 1])
        self.assertFalse(model.risk_accepts([1, 1, 1], [.99, .99, .99]))
        self.assertTrue(model.risk_accepts([0, 0, 0], [0, 0, 0]))

    def test_risk_nonfinite_or_missing_rejected(self):
        for bounds, tolerances in (([math.nan], [.1]), ([0], []), ([-.1], [.1]), ([0], [-.1])):
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                model.risk_accepts(bounds, tolerances)

    def test_initial_performance_required(self):
        result = model.run_design("B", dict(model.BASE), 5, lambda _: True, S0=(.449, .6))
        self.assertEqual(result["classes"], "GGGGG")
        self.assertFalse(result["protected_horizon"])
        self.assertFalse(result["retained_gain"])

    def test_initial_floor_equality(self):
        result = model.run_design("B", dict(model.BASE), 5, lambda _: True, S0=(.45, .6))
        self.assertTrue(result["protected_horizon"])
        self.assertTrue(result["retained_gain"])

    def test_fixed_baseline_prevents_ratchet(self):
        baseline = .65
        local = []
        protected = []
        for t in range(5):
            step = dict(SB=True, MQ=.8, SC=baseline - .009 * t,
                        SCn=baseline - .009 * (t + 1), SCmin=baseline - .009 * t)
            local.append(step["SCn"] >= step["SC"] - model.EPS_S)
            protected.append(model.review_protection_holds(step, baseline, .01))
        self.assertTrue(all(local))
        self.assertEqual(protected, [True, False, False, False, False])
        self.assertAlmostEqual(baseline - 5 * .009, .605)

    def test_fixed_baseline_review_equality_and_breach(self):
        step = dict(SB=True, MQ=.8, SC=.65, SCn=.64, SCmin=.64)
        self.assertTrue(model.review_protection_holds(step, .65, .01))
        self.assertFalse(model.review_protection_holds(dict(step, SCn=.64 - 1e-12), .65, .01))

    def test_no_hidden_terminal_index_tolerance(self):
        self.assertEqual(model.EPS_C, 0)
        result = model.run_design("B", dict(model.BASE), 5, lambda _: False)
        self.assertFalse(result["retained_gain"])

    def test_required_route_veto_all_monitoring_forms(self):
        for form in ("product", "geometric", "minimum"):
            result = model.period(model.DESIGNS["B"], dict(model.BASE, monitoring_form=form,
                                  available_routes=(False,)), [.4] * 6, [.55, .6], 1)
            self.assertEqual(result["MQ"], 0)
            self.assertEqual(result["ell"], 0)
            self.assertFalse(result["required_routes_pass"])

    def test_partial_required_group_exclusion(self):
        self.assertFalse(model.required_routes_hold((True, True), (True, False)))
        self.assertTrue(model.required_routes_hold((True, False), (True, False)))

    def test_accuracy_correct_estimand_and_unassessed(self):
        self.assertEqual(model.invalid_lesson_sensitivity(20, 0)["value"], 1)
        missing = model.invalid_lesson_sensitivity(0, 0)
        self.assertIsNone(missing["value"])
        self.assertEqual(missing["status"], "Unassessed_NoIndependentlyInvalidLessons")
        with self.assertRaises(ValueError):
            model.invalid_lesson_sensitivity(-1, 1)

    def test_zero_efficiency_and_nonpositive_qv(self):
        for eta, qv in ((0, .5), (.6, 0), (.6, -.1)):
            self.assertEqual(model.critical_learning(eta, qv, .1)["status"],
                             "NoPositiveGainAtAnyLearningLevel")

    def test_critical_threshold_strict_and_unattainable(self):
        self.assertEqual(model.critical_learning(.6, .5, 0)["value"], 0)
        self.assertEqual(model.critical_learning(.1, .2, .1)["status"], "UnattainableWithinUnitLearning")
        self.assertEqual(model.critical_learning(.5, .2, .1)["status"], "UnattainableWithinUnitLearning")

    def test_eta_q_exact_scale_nonidentifiability(self):
        self.assertLessEqual(numerical.identifiability()["maximum_trajectory_difference"], 1e-12)

    def test_full_analytic_capacity_jacobian(self):
        result = numerical.analytic_bistability()
        self.assertEqual([row["locally_attracting_capacity_subsystem"] for row in result["equilibria"]],
                         [True, False, True])
        self.assertAlmostEqual(result["fold_decay"], .16)
        self.assertAlmostEqual(result["equilibria"][1]["spectral_radius"], 1.0662741699796952)

    def test_raw_projection_diagnostics(self):
        result = model.period(model.DESIGNS["B"], dict(model.BASE, sig=100, cap=math.inf),
                              [.4] * 6, [.55, .6], 1)
        self.assertLess(result["raw_impact_proxy"][0], 0)
        self.assertEqual(result["Smin"][0], 0)
        self.assertEqual(result["path_evidence_status"], "ProxyOnly_NotObservedPathSafety")

    def test_horizon_exports_own_state_comparator(self):
        result = model.run_design("B", dict(model.BASE, routine=(.2, 0)), 5, lambda _: True)
        self.assertEqual(len(result["performance_cf"]), 6)
        self.assertGreater(result["performance_cf"][-1][0], result["performance_cf"][0][0])
        self.assertAlmostEqual(result["initial_performance_index"], .59)

    def test_schedule_errors_are_not_silent(self):
        for count, schedule in ((0, lambda _: True), (1, lambda _: 1)):
            with self.subTest(count=count), self.assertRaises(ValueError):
                model.run_design("B", dict(model.BASE), count, schedule)


if __name__ == "__main__":
    unittest.main()
