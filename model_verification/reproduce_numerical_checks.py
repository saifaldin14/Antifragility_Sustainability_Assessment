"""Reproduce safeguard, boundary, feedback, sampling-design and projection checks.

These are numerical/domain checks only; no field validation is performed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import subprocess

import numpy as np

import make_figures as model
from reproduce_scenario_analyses import read_log, require_verified_native_log


HERE = Path(__file__).resolve().parent
NATIVE = HERE / "run_verification.ps1"


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Expected nonempty records: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def native_checks(*, boundary=False):
    switch = "-BoundaryChecksOnly" if boundary else "-SafeguardChecksOnly"
    expected_count = 21 if boundary else 87
    run = subprocess.run(
        ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(NATIVE),
         switch], check=True, capture_output=True, text=True, encoding="utf-8-sig")
    rows = json.loads(run.stdout)
    if (len(rows) != expected_count or len({row["Id"] for row in rows}) != expected_count
            or any(row["Passed"] is not True for row in rows)):
        raise AssertionError(f"The complete native {switch} set must pass ({expected_count} checks).")
    if any(not math.isfinite(value) for row in rows for value in row["Values"]):
        raise AssertionError("Nonfinite native check output.")
    return rows


def validate_feedback_record(record):
    if record.get("closure") != "capacity_subsystem_fixed_forcing":
        raise ValueError("An explicitly exported fixed-forcing capacity closure is required.")
    for name in ("eta", "delta0", "chi", "lockR", "theta", "phi", "mcov", "macc",
                 "mtime", "mplur", "gamma", "rho"):
        model.unit(record[name], name)
    for name in ("h", "s0", "kappa"):
        model.finite(record[name], name)
    if record["h"] < 1 or record["s0"] <= 0 or record["kappa"] <= 0:
        raise ValueError("Invalid Hill/response parameters in feedback closure.")
    if record["eta"] + record["delta0"] + record["chi"] + record["lockR"] > 1:
        raise ValueError("Inadmissible feedback rate combination.")
    model.vector(record["q"], 6, "feedback relevance", normalized=True)
    model.vector(record["beta"], 6, "feedback index weights", normalized=True)
    if abs(math.fsum(record["beta"]) - 1) > 1e-12:
        raise ValueError("Feedback index weights must sum to one.")
    model.vector(record["external_capacity_change"], 6, "external capacity change")
    if any(record["external_capacity_change"]):
        raise ValueError("This stability experiment requires zero external capacity change.")
    model.required_routes_hold(record["required_routes"], record["available_routes"])
    if record["monitoring_form"] not in {"product", "geometric", "minimum"} or record["salience_form"] not in {"hill", "linear_capped"}:
        raise ValueError("Unknown feedback closure form.")
    size = len(record["signal"])
    if not size:
        raise ValueError("Nonempty feedback evidence space required.")
    model.vector(record["signal"], size, "fixed signal")
    model.vector(record["error"], size, "fixed error")
    if not record["AI"] or len(record["AI"]) != len(record["AT"]):
        raise ValueError("Matching nonempty feedback response matrices required.")
    for row in (*record["AI"], *record["AT"]):
        model.vector(row, size, "feedback response row")


def feedback_map(record, state):
    validate_feedback_record(record)
    model.vector(state, 6, "feedback capacity state", normalized=True)
    return _feedback_map(record, state)


def _feedback_map(record, state):
    """Equation-derived six-capacity map, with the exported fixed exogenous inputs."""
    c = np.asarray(state, dtype=float)
    if c.shape != (6,) or not np.isfinite(c).all() or np.any(c < 0) or np.any(c > 1):
        raise ValueError("A finite six-capacity state in [0,1] is required.")
    signal, error = np.asarray(record["signal"]), np.asarray(record["error"])
    evidence = signal + error
    strength = model.stable_norm(evidence)
    direction = evidence / strength if strength else np.zeros_like(evidence)
    valid = model.validity(signal, error)
    mod = lambda value: record["theta"] + (1 - record["theta"]) * value
    mq = model.monitoring_quality((record["mcov"], record["mtime"],
                                   record["mplur"] * mod(c[model.iD])), record["monitoring_form"])
    if not model.required_routes_hold(record["required_routes"], record["available_routes"]):
        mq = 0
    salience = (model.salience(strength, record["s0"], record["h"])
                if record["salience_form"] == "hill" else min(1, strength / record["s0"]))
    ell = mq * record["gamma"] * salience
    retained = record["rho"] * mod(c[model.iL]) * ell
    desired = ell * (np.asarray(record["AI"]) + c[model.iT] * np.asarray(record["AT"])) @ direction
    bound = record["kappa"] * (c[model.iF] + c[model.iO] + c[model.iR])
    utilization = min(1.0, model.stable_norm(desired) / bound) if bound else 0.0
    zeta = np.zeros(6)
    zeta[model.iO] = record["chi"] * utilization
    zeta[model.iR] = zeta[model.iO] + record["lockR"]
    decay = record["delta0"] * (1 - record["phi"] * record["rho"] * mod(c[model.iL]))
    result = np.asarray([
        model._capacity_coordinate(c[k], record["eta"], retained, record["q"][k],
                                   valid, record["macc"], decay, zeta[k])
        for k in range(6)
    ])
    if not np.isfinite(result).all() or np.any(result < 0) or np.any(result > 1):
        raise AssertionError("Equation-derived capacity map left its numerical domain.")
    return result


def jacobian(function, state, step):
    state = np.asarray(state, dtype=float)
    columns = []
    for k in range(6):
        plus, minus = state.copy(), state.copy()
        plus[k] = min(1, state[k] + step)
        minus[k] = max(0, state[k] - step)
        columns.append((function(plus) - function(minus)) / (plus[k] - minus[k]))
    return np.column_stack(columns)


def feedback_checks(text):
    records = [json.loads(line.split("|", 1)[1])
               for line in text.splitlines() if line.startswith("FEEDBACK|")]
    if not records or len({record["system_id"] for record in records}) != len(records):
        raise AssertionError("Complete, unique distinct-terminal-state exports are required.")
    summaries, states = [], []
    for record in records:
        validate_feedback_record(record)
        if len(record["initial_states"]) != 5 or len(record["terminal_states"]) != 5:
            raise AssertionError("Each sampled system requires five complete initial/terminal states.")
        function = lambda state: _feedback_map(record, state)
        for index, (initial, terminal) in enumerate(zip(record["initial_states"], record["terminal_states"]), 1):
            terminal = np.asarray(terminal)
            residual = float(np.max(np.abs(function(terminal) - terminal)))
            replayed = np.asarray(initial, dtype=float)
            for replay_steps in range(1, 200001):
                following = function(replayed)
                if np.max(np.abs(following - replayed)) <= 1e-13:
                    replayed = following
                    break
                replayed = following
            else:
                raise AssertionError(f"Initial-state replay did not converge: system {record['system_id']}, start {index}")
            replay_error = float(np.max(np.abs(replayed - terminal)))
            j1, j2 = jacobian(function, terminal, 1e-5), jacobian(function, terminal, 1e-6)
            radius = float(np.max(np.abs(np.linalg.eigvals(j2))))
            difference = float(np.max(np.abs(j2 - j1)))
            largest_basin_error = 0.0
            for k in range(6):
                for sign in (-1, 1):
                    perturbed = terminal.copy()
                    perturbed[k] = np.clip(perturbed[k] + sign * 1e-4, 0, 1)
                    for _ in range(10000):
                        following = function(perturbed)
                        if np.max(np.abs(following - perturbed)) <= 1e-13:
                            perturbed = following
                            break
                        perturbed = following
                    largest_basin_error = max(largest_basin_error, float(np.max(np.abs(perturbed - terminal))))
            if residual > 1e-8 or not radius < 1 or difference > 1e-3 or largest_basin_error > 1e-5 or replay_error > 1e-7:
                raise AssertionError(f"Local numerical feedback check failed: {record['system_id']}/{index}")
            summaries.append(dict(
                system_id=record["system_id"], initialization=index, fixed_point_residual=residual,
                spectral_radius=radius, jacobian_step_difference=difference,
                independent_initial_replay_error=replay_error, independent_initial_replay_steps=replay_steps,
                largest_perturbation_return_error=largest_basin_error,
                perturbations=12, perturbation_size=1e-4, maximum_return_steps=10000,
                scope="Local numerical capacity-subsystem evidence only; not full-state or global stability"))
            states.append(dict(system_id=record["system_id"], initialization=index,
                               initial=list(initial), terminal=terminal.tolist(), jacobian=j2.tolist(),
                               eigenvalues=[[float(value.real), float(value.imag)]
                                            for value in np.linalg.eigvals(j2)]))
    return records, summaries, states


def analytic_bistability():
    a, decay = 0.64, 0.08
    equilibria = (0.0, (1 - math.sqrt(1 - 4 * decay / a)) / 2,
                  (1 + math.sqrt(1 - 4 * decay / a)) / 2)
    rows = []
    for c in equilibria:
        diagonal = 1 - decay - a * c * c
        outer = a * c * (1 - c) * np.ones((6, 1)) @ np.array([[0, 1, 0, 1, 0, 0]])
        matrix = diagonal * np.eye(6) + outer
        rows.append(dict(capacity=c, spectral_radius=float(max(abs(np.linalg.eigvals(matrix)))),
                         locally_attracting_capacity_subsystem=bool(max(abs(np.linalg.eigvals(matrix))) < 1)))
    return dict(A=a, decay=decay, fold_decay=a / 4, equilibria=rows,
                scope="Exact invariant-diagonal equilibria and analytic six-state capacity Jacobian; not full performance-state attraction.")


def identifiability():
    first, second = dict(model.BASE), dict(model.BASE, eta=0.75, q=0.4)
    a = model.run_design("B", first, 30, lambda _: True)
    b = model.run_design("B", second, 30, lambda _: True)
    error = float(np.max(np.abs(np.asarray(a["components"]) - np.asarray(b["components"]))))
    if error > 1e-12:
        raise AssertionError("The exact eta-q scale symmetry was not reproduced.")
    return dict(first_eta_q=[0.6, 0.5], second_eta_q=[0.75, 0.4],
                common_composite=0.3, periods=30, maximum_trajectory_difference=error,
                inference="Do not fit eta and q separately without independent justified scale information; other confounding can remain.")


def sampling_stress():
    size, seed = 512, 20261024
    rng = np.random.default_rng(seed)
    names = ("delta0", "eta", "h", "theta", "phi", "q", "noise", "chi", "lockR",
             "s0", "kappa", "sig", "c0")
    bounds = ((.01, .10), (.2, .8), (1, 3), (0, 1), (0, 1), (.1, 1), (.05, .6),
              (0, .08), (0, .01), (.2, .8), (.1, .5), (.3, 1.2), (.2, .6))
    lhs = np.column_stack([(rng.permutation(size) + rng.random(size)) / size for _ in names])
    summaries, raw = [], []
    for label in ("independent_lhs", "positive_rank_q_theta", "negative_rank_q_theta", "coupled_resources"):
        design = lhs.copy()
        if label in ("positive_rank_q_theta", "negative_rank_q_theta"):
            order = np.argsort(design[:, names.index("theta")])
            values = np.sort(design[:, names.index("q")])
            design[order, names.index("q")] = values if label.startswith("positive") else values[::-1]
        counts = dict(all_reviews_gain=0, any_review_gain=0, positive_contrast=0,
                      index_increase=0, protected_horizon=0, retained_horizon_gain=0)
        for index, row in enumerate(design, 1):
            p = dict(model.BASE, **{key: lo + (hi - lo) * value for key, (lo, hi), value in zip(names, bounds, row)})
            resource = 0.25 + 0.75 * row[names.index("q")] if label == "coupled_resources" else 1.0
            governance = dict(model.DESIGNS["B"])
            governance["gamma"] *= resource
            governance["rho"] *= resource
            c, cf = [p["c0"]] * 6, [p["c0"]] * 6
            s, sf, classes = [0.55, 0.60], [0.55, 0.60], []
            protected = model.protection_holds(s, (0.45, 0), c, p["floors"], (), ())
            for _ in range(5):
                actual, counter = model.period(governance, p, c, s, 1), model.period(governance, p, c, s, 0)
                counter["C"] = c
                classes.append(model.latent_class_of(actual, counter))
                safe = model.review_protection_holds(actual, .2 * .55 + .8 * .60, p["horizon_tolerance"])
                protected = protected and safe
                c, s = actual["Cn"], actual["Sn"]
                comparator = model.period(governance, p, cf, sf, 0)
                cf, sf = comparator["Cn"], comparator["Sn"]
            final, contrast = math.fsum(c) / 6, (math.fsum(c) - math.fsum(cf)) / 6
            sequence = "".join(classes)
            retained = protected and contrast > 0 and final >= p["c0"]
            counts["all_reviews_gain"] += int(sequence == "GGGGG")
            counts["any_review_gain"] += int("G" in sequence)
            counts["positive_contrast"] += int(contrast > 0)
            counts["index_increase"] += int(final > p["c0"])
            counts["protected_horizon"] += int(protected)
            counts["retained_horizon_gain"] += int(retained)
            raw.append(dict(sampling_design=label, draw=index, **{name: p[name] for name in names},
                            resource_multiplier=resource, final_index=final, contrast=contrast,
                            classes=sequence, protected_horizon=protected, retained_horizon_gain=retained))
        summaries.append(dict(sampling_design=label, draws=size, seed=seed, **counts,
                              interpretation="Artificial design stress, not empirical probabilities, priors or prevalence"))
    return summaries, raw


def cross_language_boundary_checks(cases):
    by_id = {row["Id"]: row for row in cases}
    fused = model.canonical_evidence([.2, .6, .2], [1] * 3, [0] * 3, [1] * 3,
                                     ["a:v1", "b:v1", "a:v1"])[0]
    expected = by_id["AC01_mixed_source_copy"]["Values"][1]
    initial = model.run_design("B", dict(model.BASE), 5, lambda _: True, S0=(.449, .60))
    social_error = abs(initial["performance"][1][0] - by_id["AC14_initial_performance"]["Values"][1])
    if fused != expected or initial["protected_horizon"] or initial["retained_gain"] or social_error > 1e-12:
        raise AssertionError("Cross-language fusion/initial-window boundary disagreement.")
    return dict(mixed_source_fused=fused, native_fused=expected,
                initial_below_floor_classes=initial["classes"],
                initial_below_floor_protected=initial["protected_horizon"],
                first_review_social_cross_language_error=social_error,
                stable_large_validity=model.validity([1e308], [1e308]),
                stable_tiny_hill=model.salience(1e-200, 1e-200, 2),
                accuracy=model.invalid_lesson_sensitivity(20, 0),
                accuracy_without_invalid_lessons=model.invalid_lesson_sensitivity(0, 0),
                zero_eta_threshold=model.critical_learning(0, .5, 0),
                native_raw_projection_probe_first_component=by_id["AC18_raw_projection"]["Values"])


def projection_records():
    rows = []
    for signal in (.8, 100):
        for name in ("A", "B"):
            result = model.period(model.DESIGNS[name], dict(model.BASE, sig=signal, cap=math.inf),
                                  [.4] * 6, [.55, .6], 1)
            row = dict(design=name, disturbance=signal, safeguard_pass=result["SB"],
                       path_evidence_status=result["path_evidence_status"])
            for k in range(2):
                row[f"raw_performance_{k}"] = result["raw_performance"][k]
                row[f"projected_performance_{k}"] = result["Sn"][k]
                row[f"raw_impact_proxy_{k}"] = result["raw_impact_proxy"][k]
                row[f"projected_impact_proxy_{k}"] = result["Smin"][k]
            for k, label in enumerate(("F", "D", "O", "L", "R", "T")):
                row[f"raw_{label}"] = result["raw_capacity"][k]
                row[f"projected_{label}"] = result["Cn"][k]
            rows.append(row)
    return rows


def reproduce(log, output):
    text = read_log(log)
    require_verified_native_log(text)
    output.mkdir(parents=True, exist_ok=True)
    cases = native_checks()
    boundary_cases = native_checks(boundary=True)
    records, checks, states = feedback_checks(text)
    summaries, draws = sampling_stress()
    write_json(output / "Safeguard_Checks.json", cases)
    write_json(output / "Boundary_Checks.json", boundary_cases)
    write_json(output / "Feedback_Sampled_Maps.json", records)
    write_csv(output / "Feedback_Local_Checks.csv", checks)
    write_json(output / "Feedback_Jacobians_and_States.json", states)
    write_csv(output / "Sampling_Design_Stress.csv", summaries)
    write_csv(output / "Sampling_Design_Draws.csv", draws)
    write_csv(output / "Projection_Diagnostics.csv", projection_records())
    result = dict(
        status="PASS", native_safeguard_checks=len(cases), native_boundary_checks=len(boundary_cases),
        distinct_terminal_state_systems=len(records), terminal_states_checked=len(checks),
        capacity_perturbation_checks=sum(row["perturbations"] for row in checks),
        maximum_fixed_point_residual=max(row["fixed_point_residual"] for row in checks),
        maximum_capacity_spectral_radius=max(row["spectral_radius"] for row in checks),
        maximum_independent_initial_replay_error=max(row["independent_initial_replay_error"] for row in checks),
        maximum_local_perturbation_return_error=max(row["largest_perturbation_return_error"] for row in checks),
        identifiability=identifiability(), analytic_bistability=analytic_bistability(),
        cross_language_boundary_checks=cross_language_boundary_checks(cases), projection_diagnostic_rows=4,
        sampling_stress=summaries, artificial_design_runs=len(draws),
        zero_of_59_independent_bernoulli_failures_one_sided_95_upper=1 - 0.05 ** (1 / 59),
        evidence_status=model.SIMULATION_EVIDENCE_STATUS,
        interpretation="Numerical/domain/theory checks only; no measured capacity validation, hazard-coverage validation or independent human authorization.")
    write_json(output / "Numerical_Checks_Summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-log", type=Path, default=HERE / "Verified_Native_Output.txt")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    print(json.dumps(reproduce(args.input_log, args.output), indent=2, allow_nan=False))
