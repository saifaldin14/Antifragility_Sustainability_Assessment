"""Reproduce documentary sufficiency, rating provenance, and agreement statistics."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np


CASES = {
    "C1": "evolv1", "C2": "Ahrens Street West",
    "C3": "One Hundred Victoria Street South",
    "C4": "Northfield Station Area", "C5": "Breithaupt Block",
}
CRITERIA = {
    "G1": "Dated decision unit, prior reference and material disturbance distinct from routine learning",
    "G2": "Observed response and evidence-to-decision mechanism",
    "G3": "Comparable baseline/follow-up measurements of six capacities",
    "G4": "Retention or transfer beyond immediate correction",
    "G5": "Credible counterfactual and rival explanations",
    "G6": "Independent performance, group, temporary, cumulative and critical-capacity protections over a declared horizon",
    "G7": "Capability-contrast measurement and evidence sufficiency",
    "G8": "Source period, phase, and applicability established",
}
STATUSES = {"E", "P", "U", "X", "N"}
DECLARED_REDUCTIONS = {"C2-D02", "C4-D02", "C5-D02", "C5-D12", "C5-D14", "C5-D15"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Missing or duplicate CSV headers: {path.name}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty CSV: {path.name}")
    for number, row in enumerate(rows, 2):
        if set(row) != set(fields) or any(value is None for value in row.values()):
            raise ValueError(f"CSV field-count error: {path.name}:{number}")
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty result table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def ratings(values) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not len(array) or not np.all(np.isfinite(array)):
        raise ValueError("Ratings must be a nonempty finite one-dimensional array")
    if np.any(array != np.floor(array)) or np.any(array < 1) or np.any(array > 6):
        raise ValueError("Ratings must be integers 1 through 6")
    return array.astype(np.int64)


def contingency(first, second) -> np.ndarray:
    a, b = ratings(first), ratings(second)
    if a.shape != b.shape:
        raise ValueError("Paired ratings must have equal lengths")
    table = np.zeros((6, 6))
    np.add.at(table, (a - 1, b - 1), 1)
    return table


def weighted_kappa(first, second, power: int = 2) -> float | None:
    if power not in (1, 2):
        raise ValueError("Only linear and quadratic weighting are supported")
    observed = contingency(first, second)
    observed /= observed.sum()
    expected = np.outer(observed.sum(axis=1), observed.sum(axis=0))
    levels = np.arange(6)
    distance = np.abs(levels[:, None] - levels[None, :]) ** power
    denominator = float(np.sum(expected * distance))
    return 1 - float(np.sum(observed * distance)) / denominator if denominator > 0 else None


def ordinal_alpha(first, second) -> float | None:
    table = contingency(first, second)
    coincidence = table + table.T
    marginal = coincidence.sum(axis=1)
    midranks = np.cumsum(marginal) - marginal / 2
    distance = (midranks[:, None] - midranks[None, :]) ** 2
    expected = float(np.sum(np.outer(marginal, marginal) * distance))
    observed = float(np.sum(coincidence * distance))
    return 1 - (coincidence.sum() - 1) * observed / expected if expected > 0 else None


def bootstrap(first, second, metric, blocks=None, exact=False) -> dict:
    a, b = ratings(first), ratings(second)
    contingency(a, b)
    if exact and blocks is None:
        raise ValueError("Exact enumeration requires case blocks")
    groups = blocks if blocks is not None else [np.array([i]) for i in range(len(a))]
    if not groups or any(len(group) == 0 for group in groups):
        raise ValueError("Bootstrap groups must be nonempty")
    flat = np.concatenate(groups)
    if sorted(flat.tolist()) != list(range(len(a))):
        raise ValueError("Bootstrap groups must partition the paired observations")
    rng = np.random.default_rng(20261019)
    selections = (
        itertools.product(range(len(groups)), repeat=len(groups)) if exact
        else (rng.integers(0, len(groups), len(groups)) for _ in range(2000))
    )
    values, undefined = [], 0
    for selected in selections:
        indices = np.concatenate([groups[int(i)] for i in selected])
        value = metric(a[indices], b[indices])
        if value is None:
            undefined += 1
        else:
            values.append(value)
    if not values:
        raise ValueError("Statistic undefined in every bootstrap resample")
    return {
        "percentile_95": np.percentile(values, [2.5, 97.5]).tolist(),
        "defined_resamples": len(values), "undefined_resamples": undefined,
        "method": "exact ordered case resamples" if exact else "2000 seeded percentile resamples",
        "resampling_unit": "case" if blocks is not None else "cell",
        "seed": None if exact else 20261019,
    }


def load_ratings(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    expected = {f"{case}-D{d:02}" for case in CASES for d in range(1, 16)}
    identifiers = [row["item_id"] for row in rows]
    if len(rows) != 75 or set(identifiers) != expected:
        raise ValueError("Expected exactly 75 unique case-dimension cells")
    for row in rows:
        if row["item_id"] != f"{row['case_id']}-{row['dimension_id']}":
            raise ValueError(f"Inconsistent rating identifier: {row['item_id']}")
        if row["consensus"].strip():
            raise ValueError("This unreconciled dataset has no recorded consensus ratings")
        if not row["dimension"].strip() or not row["revision_provenance"].strip():
            raise ValueError("Dimension and revision provenance are required")
    rows.sort(key=lambda row: row["item_id"])
    a, initial, final = [
        ratings([row[field] for row in rows])
        for field in ("rater1", "second_round1", "second_final")
    ]
    blocks = [np.array([i for i, row in enumerate(rows) if row["case_id"] == case]) for case in CASES]
    if any(len(block) != 15 for block in blocks):
        raise ValueError("Expected fifteen dimensions in each of five cases")
    return a, initial, final, blocks


def assess_evidence(rows: list[dict[str, str]], sources: list[dict[str, str]]) -> dict:
    source_cases = {row["source_id"]: row["case_id"] for row in sources}
    if len(source_cases) != len(sources):
        raise ValueError("Duplicate source identifier")
    cells = {}
    for row in rows:
        case, criterion, status = row["case_id"], row["criterion"], row["status"]
        if case not in CASES or criterion not in CRITERIA or status not in STATUSES:
            raise ValueError(f"Invalid case/criterion/status: {case}/{criterion}/{status}")
        key = (case, criterion)
        if key in cells:
            raise ValueError(f"Duplicate empirical-sufficiency cell: {key}")
        if not all(row[field].strip() for field in ("evidence_summary", "locator", "limitation")):
            raise ValueError(f"Missing evidence, locator, or limitation: {key}")
        references = row["source_ids"].split(";")
        if any(ref not in source_cases or source_cases[ref] not in {case, "ALL"} for ref in references):
            raise ValueError(f"Unknown or cross-case source: {key}")
        cells[key] = status
    if set(cells) != {(case, criterion) for case in CASES for criterion in CRITERIA}:
        raise ValueError("Expected exactly forty empirical-sufficiency cells")
    result = {}
    for case, name in CASES.items():
        unresolved = [criterion for criterion in CRITERIA if cells[case, criterion] != "E"]
        result[case] = {
            "case": name,
            "criteria": {criterion: cells[case, criterion] for criterion in CRITERIA},
            "criteria_not_fully_established": unresolved,
            "inference_status": (
                "InsufficientForCausalTest" if unresolved
                else "ReadyForIndependentEffectAndSafeguardAnalysis"
            ),
        }
    return result


def compare_baseline(rows, baseline) -> list[dict]:
    if len(baseline) != 15 or {row["dimension_id"] for row in baseline} != {f"D{i:02}" for i in range(1, 16)}:
        raise ValueError("Expected fifteen unique published baseline dimensions")
    lookup = {row["dimension_id"]: row for row in baseline}
    comparison = []
    for row in rows:
        original = int(lookup[row["dimension_id"]][row["case_id"]])
        ratings([original])
        change = int(row["rater1"]) - original
        comparison.append({
            "item_id": row["item_id"], "published_baseline": original,
            "current_rater1": int(row["rater1"]), "change": change,
            "declared_reduction": row["item_id"] in DECLARED_REDUCTIONS,
            "source_id": "T1", "source_locator": lookup[row["dimension_id"]]["locator"],
        })
    changed = {row["item_id"] for row in comparison if row["change"]}
    if changed != DECLARED_REDUCTIONS or any(row["change"] not in (-1, 0) for row in comparison):
        raise ValueError("Baseline changes differ from the six declared one-level reductions")
    return comparison


def describe(first, second, blocks) -> dict:
    table = contingency(first, second) / len(first)
    expected = np.outer(table.sum(axis=1), table.sum(axis=0))
    nearby = np.abs(np.arange(6)[:, None] - np.arange(6)[None, :]) <= 1
    return {
        "ordinal_alpha": ordinal_alpha(first, second),
        "quadratic_kappa": weighted_kappa(first, second),
        "linear_kappa": weighted_kappa(first, second, 1),
        "exact_agreements": int(np.count_nonzero(first == second)),
        "within_one_category": int(np.count_nonzero(np.abs(first - second) <= 1)),
        "mean_second_minus_first": float(np.mean(second - first)),
        "marginal_chance_exact": float(np.trace(expected)),
        "marginal_chance_within_one": float(expected[nearby].sum()),
        "cell_bootstrap_alpha": bootstrap(first, second, ordinal_alpha),
        "cell_bootstrap_quadratic_kappa": bootstrap(first, second, weighted_kappa),
        "seeded_case_bootstrap_alpha": bootstrap(first, second, ordinal_alpha, blocks),
        "exact_case_bootstrap_alpha": bootstrap(first, second, ordinal_alpha, blocks, exact=True),
        "by_case": {
            case: {
                "ordinal_alpha": ordinal_alpha(first[block], second[block]),
                "quadratic_kappa": weighted_kappa(first[block], second[block]),
                "exact_agreements": int(np.count_nonzero(first[block] == second[block])),
                "second_rater_distinct_levels": len(np.unique(second[block])),
            } for case, block in zip(CASES, blocks)
        },
        "by_dimension": {
            f"D{d:02}": {
                "ordinal_alpha": ordinal_alpha(first[d - 1::15], second[d - 1::15]),
                "quadratic_kappa": weighted_kappa(first[d - 1::15], second[d - 1::15]),
            } for d in range(1, 16)
        },
    }


def reproduce(input_dir: Path, output_dir: Path) -> dict:
    rows = read_csv(input_dir / "Ratings.csv")
    first, initial, final, blocks = load_ratings(rows)
    comparison = compare_baseline(rows, read_csv(input_dir / "Thesis_Table20_Ratings.csv"))
    evidence = assess_evidence(read_csv(input_dir / "Case_Criteria.csv"),
                               read_csv(input_dir / "Source_Ledger.csv"))
    changed = final != initial
    result = {
        "paired_cells": len(first), "case_clusters": len(blocks),
        "first_round": describe(first, initial, blocks),
        "final_round": describe(first, final, blocks),
        "second_rater_revisions": int(np.count_nonzero(changed)),
        "revisions_upward": int(np.count_nonzero(final > initial)),
        "revisions_toward_first": int(np.count_nonzero(changed & (np.abs(final - first) < np.abs(initial - first)))),
        "recorded_consensus_cells": 0,
        "primary_uncertainty_method": "Exact enumeration of all 3,125 ordered five-case block resamples; final ordinal-alpha percentile interval.",
        "sensitivity_uncertainty_methods": ["Cell bootstrap", "Seeded 2,000-resample case bootstrap"],
        "uncertainty_limit": "Descriptive with five related purposive case clusters; enumeration does not establish nominal coverage or blinded reliability.",
        "provenance": "Nonblind ratings; final second ratings revised after agreement results; rationales added afterwards.",
        "interpretation": "Diagnostic integration judgments, not six-capacity observations or causal validation.",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, content in (("Rating_Statistics.json", result), ("Case_Empirical_Status.json", evidence)):
        (output_dir / name).write_text(json.dumps(content, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(output_dir / "Published_Rating_Comparison.csv", comparison)
    write_csv(output_dir / "Second_Rater_Revision_History.csv", [{
        "item_id": row["item_id"], "rater1": int(a),
        "second_round1": int(b), "second_final": int(c), "change": int(c - b),
        "consensus": "", "revision_provenance": row["revision_provenance"],
    } for row, a, b, c in zip(rows, first, initial, final)])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parent / "data")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    results = reproduce(args.input, args.output)
    print(json.dumps({
        "paired_cells": results["paired_cells"],
        "final_alpha": results["final_round"]["ordinal_alpha"],
        "exact_case_interval": results["final_round"]["exact_case_bootstrap_alpha"]["percentile_95"],
        "second_rater_revisions": results["second_rater_revisions"],
        "output": args.output.name,
    }, indent=2))
