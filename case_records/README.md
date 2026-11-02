# Case records (Supplementary Material 2)

Documentary assessment of five Waterloo Region cases: integration ratings on
fifteen sustainability dimensions by two raters, and evidence-sufficiency
coding on eight criteria (G1-G8). The cases are C1 evolv1, C2 Ahrens Street
West, C3 One Hundred Victoria Street South, C4 Northfield Station Area and C5
Breithaupt Block. The assessment uses published documents only. No field
measurements or interviews were conducted, and it is not a causal validation of
the dynamic model.

## Running

Tested with Python 3.12.10 and NumPy 2.5.4 (`requirements.txt`). From this
folder:

```powershell
python -m pip install -r requirements.txt pytest
python -m pytest -q -p no:cacheprovider        # 24 passed
python -B reproduce_cases.py --output .\recomputed
```

`python -B reproduce_cases.py` without `--output` rewrites `results/`. The
tests can also be run with `python -B -m unittest -v test_reproduction.py`.
Regenerated files are byte-identical to those in `results/`.

## Files and schemas

### `data/` (inputs)

- `Ratings.csv` - 75 cells (5 cases x 15 dimensions D01-D15). `rater1`: current
  Rater 1 score; `second_round1`: initial second rating; `second_final`: final
  second rating; `second_intermediate` and `intermediate_note`: a documented
  intermediate second rating; `second_confidence`, `second_locator`,
  `second_rationale`: final evidence locator and rationale; `consensus`: always
  empty (no consensus was recorded); `revision_provenance`: how the final second
  rating was reached. Scores are integers on the six-level integration rubric:
  1 None, 2 Symbolic, 3 Limited, 4 Moderate, 5 Considerable, 6 Full.
- `Thesis_Table20_Ratings.csv` - the published five-case matrix (thesis Table 20,
  source T1). Six current Rater 1 cells are one level lower; 69 are unchanged. A
  prose example in the thesis gives evolv1 resilience as 6, but the published
  table gives 4, which matches the current rating.
- `Case_Criteria.csv` - 40 cells (5 cases x G1-G8): status code, evidence
  summary, source IDs, locators and limitation.
- `Source_Ledger.csv` - the 13 case-related sources (including the thesis) and
  two contextual municipal guidance documents (`source_id` C1/C2; these are not
  case outcome evidence). `source_id` is distinct from `case_id`. Gives URLs,
  source period, locators, access status and limitations.
- `Rationale_Qualification_Gaps.csv` - four cells (C1-D05, C1-D09, C2-D04,
  C2-D05) whose rationale does not by itself justify the selected level. Scores
  and rationales are left unchanged; only an independent source/rubric
  assessment could settle the level.
- `Source_Archive_Metadata.csv` - SHA-256 fingerprints, sizes and page counts of
  two downloaded sources. The documents are not redistributed.

### `results/` (produced by `reproduce_cases.py`)

- `Rating_Statistics.json` - agreement between Rater 1 and the initial and final
  second ratings: ordinal Krippendorff's alpha, quadratic and linear weighted
  kappa, exact and within-one agreement, mean difference, chance agreement,
  bootstrap and exact-enumeration intervals, and by-case and by-dimension values.
  Undefined coefficients are JSON `null`; undefined resamples are counted.
- `Case_Empirical_Status.json` - G1-G8 codes and the resulting sufficiency
  status for each case.
- `Published_Rating_Comparison.csv` - all 75 differences between the published
  matrix and current Rater 1 scores.
- `Second_Rater_Revision_History.csv` - initial and final second ratings for all
  75 cells, including unchanged cells.

## Interpretation

**Ratings.** Ratings measure documented integration of commitments, not
sustainability performance or the model's six capacities. Rater 2 (a member of
the author team) rated nonblind. Final second ratings were revised after the
agreement results had been disclosed (30 cells changed: 29 upward, 21 toward
Rater 1), and final locators and rationales were recorded afterwards. No
consensus or blinded independent assessment is claimed.

**Uncertainty.** The primary uncertainty description is exact enumeration of
all 3,125 ordered five-case block resamples: final ordinal alpha 0.444, with a
percentile interval from -0.216 to 0.680. With five related, purposively chosen
cases this is descriptive and does not establish nominal coverage. The cell
bootstrap and the seeded case bootstrap (2,000 resamples, seed 20261019) are
sensitivity outputs only. Enumeration removes Monte Carlo error, not the
inferential limits of five cases.

**Evidence-sufficiency codes.** E = established for the stated criterion;
P = partial; U = unestablished in this corpus; X = contradictory; N = not
applicable. U does not mean zero effect, zero capability or project failure.
G1-G8 cover: a dated decision unit, prior reference and material disturbance
distinct from routine learning (G1 requires a predeclared materiality rule;
dated discrepancies alone do not meet it); observed response and mechanism;
comparable measurements of the six capacities; retention or transfer;
counterfactual and rival explanations; independent performance, group,
temporary, cumulative and critical-capacity protections over a declared
horizon; measurement sufficiency; and source period/phase.

**What the statuses mean.** Every case is `InsufficientForCausalTest`: none has
enough documentary evidence for a full test of attributable capability gain
under the safeguards. Even eight E codes would only give
`ReadyForIndependentEffectAndSafeguardAnalysis`, never a verified
antifragility label or a positive effect. The coding is nonblind and not
preregistered; AI tool use is disclosed in the article.

**Sources.** The thesis (T1: Mohamed, 2026, University of Waterloo,
https://hdl.handle.net/10012/24396) shares research provenance with this work
and is not independent replication; its 707-page PDF was indexed and selected
passages and Table 20 were examined. The thesis phrase "stakeholder interviews"
(Section 4.4.2, printed p. 85) is not treated as interview data. Repeated
narratives and syndicated reports can share provenance, so multiple URLs are not
independent corroboration. For Breithaupt Block, the current B4 page calls
Phase III six storeys while the dated B2 release describes it as eleven; no
as-built phase crosswalk was established. Reported outreach counts are reach,
not demonstrated benefits.

**Privacy and copyright.** No copyrighted source PDFs, resident contact details
or participant photographs are included; primary-source URLs, locators and
short analytical paraphrases are provided instead.
