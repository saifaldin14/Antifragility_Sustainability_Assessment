"""Reproduce the global-sensitivity, structural, matched-benchmark, dossier and scenario analyses.

All analyses are synthetic; no field calibration is performed.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, ROUND_HALF_UP
from fractions import Fraction
import io
import hashlib
import json
import math
from pathlib import Path
import subprocess

import make_figures as model


HERE = Path(__file__).resolve().parent
NATIVE = HERE / "run_verification.ps1"


def read_log(path):
    content = path.read_bytes()
    encoding = "utf-16" if content.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8-sig"
    return content.decode(encoding)


def require_native_provenance(text):
    hashes = [line.split("|", 1)[1] for line in text.splitlines() if line.startswith("NATIVE_SOURCE_SHA256|")]
    if hashes != [hashlib.sha256(NATIVE.read_bytes()).hexdigest()]:
        raise ValueError("Native log/source SHA-256 mismatch or missing provenance; rerun full native verification.")


def require_verified_native_log(text):
    require_native_provenance(text)
    required = ("Reported values reproduced: all 64", "Test case counts as reported: all 39",
                "Scenario checks passed: 31", "Safeguard checks passed: 87",
                "Boundary checks passed: 21", "RESULT: PASS")
    if any(item not in text for item in required) or "RESULT: FAIL" in text:
        raise ValueError("A complete passing current native log is required.")


def extract_draws(text):
    lines = [line[3:] for line in text.splitlines() if line.startswith("GD|")]
    rows = list(csv.DictReader(io.StringIO("\n".join(lines))))
    if len(rows) != 2000 or [int(row["draw"]) for row in rows] != list(range(1, 2001)):
        raise ValueError("Exactly the 2,000 uniquely ordered native sensitivity draws are required.")
    expected = {"draw", "delta0", "eta", "h", "theta", "phi", "q", "noise", "chi", "lock",
                "sigma0", "kappa", "disturbance", "c0", "A_capacity5", "B_capacity5",
                "B_counterfactual5", "B_counterfactual_no_error5", "A_classes", "B_classes",
                "horizon_tolerance", "A_protected_horizon", "B_protected_horizon",
                "A_retained_gain", "B_retained_gain"}
    if any(set(row) != expected or any(value is None or value == "" for value in row.values()) for row in rows):
        raise ValueError("Every sensitivity draw requires the complete 25-column export schema.")
    return rows


def check_global(rows):
    if len(rows) != 2000:
        raise ValueError("The complete 2,000-draw experiment is required.")
    aliases = {"lock": "lockR", "sigma0": "s0", "disturbance": "sig"}
    parameters = ("delta0", "eta", "h", "theta", "phi", "q", "noise", "chi", "lock",
                  "sigma0", "kappa", "disturbance", "c0")
    counts = dict(BaboveA=0, BaboveCounterfactual=0, BcapabilityGainEveryYear=0,
                  BcapabilityGainAnyYear=0, BlevelRises=0, AanyFavourable=0,
                  BaboveNoiseFreeCounterfactual=0, BprotectedHorizon=0, BretainedHorizonGain=0)
    maximum_error = 0.0
    for row in rows:
        p = dict(model.BASE)
        p["horizon_tolerance"] = float(row["horizon_tolerance"])
        for name in parameters:
            value = float(row[name])
            if not math.isfinite(value):
                raise ValueError(f"Nonfinite parameter in draw {row['draw']}: {name}")
            p[aliases.get(name, name)] = value
        a = model.run_design("A", p, 5, lambda _: True)
        b = model.run_design("B", p, 5, lambda _: True)
        b0 = model.run_design("B", dict(p, noise=0.0), 5, lambda _: False)
        comparisons = ((a["idx"][-1], "A_capacity5"), (b["idx"][-1], "B_capacity5"),
                       (b["idxcf"][-1], "B_counterfactual5"),
                       (b0["idx"][-1], "B_counterfactual_no_error5"))
        for calculated, column in comparisons:
            expected = float(row[column])
            error = abs(calculated - expected)
            maximum_error = max(maximum_error, error)
            if not math.isfinite(expected) or error > 1e-12:
                raise AssertionError(f"Cross-language mismatch: draw {row['draw']}, {column}: {calculated}, {expected}")
        if a["classes"] != row["A_classes"] or b["classes"] != row["B_classes"]:
            raise AssertionError(f"Cross-language class mismatch in draw {row['draw']}")
        for result, prefix in ((a, "A"), (b, "B")):
            for key, column in (("protected_horizon", "protected_horizon"), ("retained_gain", "retained_gain")):
                supplied = row[f"{prefix}_{column}"]
                if supplied not in {"True", "False"}:
                    raise ValueError("Horizon exports must be explicit True/False fields.")
                if result[key] != (supplied == "True"):
                    raise AssertionError(f"Cross-language horizon mismatch: draw {row['draw']}/{prefix}/{key}")
        counts["BaboveA"] += b["idx"][-1] > a["idx"][-1]
        counts["BaboveCounterfactual"] += b["idx"][-1] > b["idxcf"][-1]
        counts["BcapabilityGainEveryYear"] += b["classes"] == "GGGGG"
        counts["BcapabilityGainAnyYear"] += "G" in b["classes"]
        counts["BlevelRises"] += b["idx"][-1] > b["idx"][0]
        counts["AanyFavourable"] += any(letter in "GRS" for letter in a["classes"])
        counts["BaboveNoiseFreeCounterfactual"] += b["idx"][-1] > b0["idx"][-1]
        counts["BprotectedHorizon"] += b["protected_horizon"]
        counts["BretainedHorizonGain"] += b["retained_gain"]
    return dict(status="PASS", draws=len(rows), numeric_comparisons=4 * len(rows),
                class_sequence_comparisons=2 * len(rows), absolute_tolerance=1e-12,
                horizon_boolean_comparisons=4 * len(rows), exported_columns=25,
                maximum_absolute_difference=maximum_error, counts=counts,
                percentages={name: float((Decimal(count) * 100 / Decimal(len(rows))).quantize(
                    Decimal("0.1"), rounding=ROUND_HALF_UP)) for name, count in counts.items()},
                interpretation="Cross-language numerical reproduction, not independent human validation or empirical prevalence.")


def native_checks():
    run = subprocess.run(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                          str(NATIVE), "-ScenarioChecksOnly"], check=True, capture_output=True,
                         text=True, encoding="utf-8-sig")
    rows = json.loads(run.stdout)
    expected_ids = ({f"SC{i:02}" for i in range(1, 20)}
                    | {f"SC20_{form}" for form in ("product", "geometric", "minimum")}
                    | {"SC21_hill", "SC21_linear_capped"}
                    | {f"SC{i:02}" for i in range(22, 29)})
    if len(rows) != len(expected_ids) or {row["Id"] for row in rows} != expected_ids:
        raise AssertionError("Incomplete or duplicated scenario checks.")
    if any(row["Passed"] is not True for row in rows):
        raise AssertionError("A native scenario check failed.")
    independent = {
        "SC05": [Fraction(3, 4), 0, 1],
        "SC08": [0, Fraction(1, 4), Fraction(1, 4)],
        "SC09": [Fraction(1, 8), Fraction(1, 8), Fraction(1, 4)],
        "SC13": [Fraction(1, 2), Fraction(11, 20), Fraction(1, 2), Fraction(3, 10), Fraction(2, 5)],
        "SC16": [Fraction(1, 5)] * 3,
        "SC17": [Fraction(1, 5), Fraction(1, 3)],
        "SC24": [Fraction(-12, 125), Fraction(18, 125)],
        "SC27": [Fraction(6, 5), 1],
        "SC28": [Fraction(9, 10), Fraction(4, 5)],
    }
    by_id = {row["Id"]: row for row in rows}
    for case_id, expected in independent.items():
        actual = by_id[case_id]["Values"]
        if len(actual) != len(expected) or any(abs(a - float(b)) > 1e-12 for a, b in zip(actual, expected)):
            raise AssertionError(f"Exact-arithmetic comparison failed: {case_id}")
    return rows


def structural_rows():
    variants = [
        ("zero_routine_boundary", {}, 1.0),
        ("routine_0.2", {"routine": (0.2, 0.0)}, 1.0),
        ("routine_0.4", {"routine": (0.4, 0.0)}, 1.0),
        ("routine_0.8", {"routine": (0.8, 0.0)}, 1.0),
        ("routine_no_error", {"routine": (0.2, 0.0), "noise": 0.0}, 1.0),
        ("geometric_monitoring", {"routine": (0.2, 0.0), "monitoring_form": "geometric"}, 1.0),
        ("minimum_monitoring", {"routine": (0.2, 0.0), "monitoring_form": "minimum"}, 1.0),
        ("concave_hill_h1", {"routine": (0.2, 0.0), "h": 1.0}, 1.0),
        ("capped_linear_salience", {"routine": (0.2, 0.0), "salience_form": "linear_capped"}, 1.0),
        ("coupled_resource_half", {"routine": (0.2, 0.0)}, 0.5),
        ("coupled_resource_quarter", {"routine": (0.2, 0.0)}, 0.25),
        ("R_floor_0.2", {"routine": (0.2, 0.0), "floors": (0, 0, 0, 0, 0.2, 0)}, 1.0),
        ("R_floor_0.4", {"routine": (0.2, 0.0), "floors": (0, 0, 0, 0, 0.4, 0)}, 1.0),
    ]
    rows = []
    for name, changes, resource in variants:
        p = dict(model.BASE, **changes)
        design = dict(model.DESIGNS["B"])
        design["gamma"] *= resource
        design["rho"] *= resource
        c = [p["c0"]] * model.K
        cf = list(c)
        s = [0.55, 0.60]
        classes = []
        for _ in range(5):
            step = model.period(design, p, c, s, 1)
            counter = model.period(design, p, c, s, 0)
            counter["C"] = c
            classes.append(model.latent_class_of(step, counter))
            c, s = step["Cn"], step["Sn"]
            cf = model.period(design, p, cf, [0.55, 0.60], 0)["Cn"]
        row = dict(scenario=name, routine_input=p["routine"][0], monitoring_form=p["monitoring_form"],
                   salience_form=p["salience_form"], h=p["h"], error=p["noise"],
                   coupled_resource_multiplier=resource, reversibility_floor=p["floors"][model.iR],
                   capacity5=sum(c) / model.K, comparator5=sum(cf) / model.K,
                   contrast5=(sum(c) - sum(cf)) / model.K, classes="".join(classes),
                   evidence_status=model.SIMULATION_EVIDENCE_STATUS)
        row.update(zip(("F", "D", "O", "L", "R", "T"), c))
        rows.append(row)
    return rows


def matched_benchmarks():
    specs = (
        ("routine_only", 0.4, [0.5] * 6, 0.5, (0.2, 0.2), (1.0, 1.0), 0.30, "Robust"),
        ("avoided_loss", 0.8, [0.8] * 6, 0.7, (0.2, 0.2), (1.0, 1.0), 0.30, "CapabilityGain"),
        ("lost_retention", 0.4, [0.35] * 6, 0.30, (0.2, 0.2), (1.0, 1.0), 0.30, "Robust"),
        ("group_breach_control", 0.4, [0.5] * 6, 0.4, (1.5, 0.0), (1.0, 1.0), 0.30, "Unacceptable"),
        ("critical_capacity_control", 0.5, [0.6, 0.6, 0.6, 0.6, 0.3, 0.6], 0.5,
         (0.2, 0.2), (1.0, 1.0), 0.40, "Unacceptable"),
    )
    rows = []
    for name, initial, final, counter, harm, limits, r_floor, expected in specs:
        floors = (0, 0, 0, 0, r_floor, 0)
        safe = model.protection_holds((0.6,), (0.4,), final, floors, harm, limits)
        response = dict(Cn=final, SB=safe, MQ=1.0, SC=0.6, SCn=0.6, SCmin=0.6,
                        d=0.0 if name == "routine_only" else 0.2)
        zero = dict(C=[initial] * 6, Cn=[counter] * 6)
        verdict = model.LATENT_CLASS_NAMES[model.latent_class_of(response, zero)]
        if verdict != expected:
            raise AssertionError(f"Constructed benchmark failed its explicit requirement: {name}")
        initial_index = sum([initial] * 6) / 6
        comparator_index = sum([counter] * 6) / 6
        rows.append(dict(scenario=name, initial_index=initial_index, final_index=sum(final) / 6,
                         comparator_index=comparator_index, contrast=sum(final) / 6 - comparator_index,
                         index_change=sum(final) / 6 - initial_index, safeguards=safe,
                         verdict=verdict, resource_units=1.0, information="same records for both assessment approaches",
                         source_type="Constructed assessment inputs; not field observations or fitted trajectories",
                         evidence_status=model.SIMULATION_EVIDENCE_STATUS))
    return rows


def dossier_rows():
    p = dict(model.BASE, routine=(0.2, 0.0))
    result = model.run_design("B", p, 5, lambda _: True)
    c, cf, s = [0.4] * 6, [0.4] * 6, [0.55, 0.60]
    burden = 0.0
    rows = []
    for year in range(1, 6):
        step = model.period(model.DESIGNS["B"], p, c, s, 1)
        cf = model.period(model.DESIGNS["B"], p, cf, [0.55, 0.60], 0)["Cn"]
        # This separately declared additive ledger is not inferred from endpoint rates.
        burden += 0.24
        safe = model.protection_holds(step["Sn"], (0.45, 0), step["Cn"], p["floors"],
                                      (0.24,), (0.30,), cumulative=(burden,), cumulative_limits=(1.0,))
        rows.append(dict(scenario="Fictional_Programme_P", forecast_review=year,
                         capacity_index=sum(step["Cn"]) / 6, comparator_index=sum(cf) / 6,
                         reversibility=step["Cn"][model.iR], reversibility_floor=0.30,
                         integrated_burden_increment=0.24, cumulative_burden=burden, cumulative_budget=1.0,
                         all_declared_protections=safe,
                         instruction=("Reject/redesign before the forecast breach" if not safe else
                                      "Escalate before next planned increment" if year == 4 else "Conditional forecast only"),
                         source_type="Constructed forecast; no observed exposure or authorized real experiment",
                         evidence_status=model.SIMULATION_EVIDENCE_STATUS))
        c, s = step["Cn"], step["Sn"]
    if rows[-1]["all_declared_protections"] or not math.isclose(burden, 1.2, abs_tol=1e-12):
        raise AssertionError("Cumulative breach must override the positive capability contrast.")
    if abs(rows[-1]["capacity_index"] - result["idx"][-1]) > 1e-12:
        raise AssertionError("Dossier and worked-example model diverged.")
    return rows


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Cannot export empty analysis: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-log", type=Path, default=HERE / "Verified_Native_Output.txt")
    parser.add_argument("--output", type=Path, default=Path("recomputed"))
    args = parser.parse_args()
    text = read_log(args.input_log)
    require_verified_native_log(text)
    draws = extract_draws(text)
    cross = check_global(draws)
    cases = native_checks()
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "Global_Sensitivity_Draws.csv", draws)
    write_csv(args.output / "Structural_Robustness.csv", structural_rows())
    write_csv(args.output / "Matched_Decision_Benchmarks.csv", matched_benchmarks())
    write_csv(args.output / "Illustrative_Dossier.csv", dossier_rows())
    for name, data in (("Model_Cross_Implementation_Checks.json", cross), ("Scenario_Checks.json", cases)):
        (args.output / name).write_text(json.dumps(data, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(cross_implementation=cross, scenario_checks=len(cases),
                          structural_variants=13, matched_benchmarks=5, dossier_forecasts=5), indent=2))


if __name__ == "__main__":
    main()
