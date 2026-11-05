# Model and verification (Supplementary Material 1)

Code, tests and outputs for the illustrative dynamic learning-capacity model:
six capacities (flexibility F, diversity D, optionality O, learning capacity L,
reversibility R, transformative capacity T), a routine-aware counterfactual,
retention, and non-compensable safeguards. The model is illustrative and
uncalibrated; every simulated result has evidence status
`NotAssessed_IllustrativeSimulation`. A simulated gain, a positive contrast, or
a complete evidence checklist is not an empirical antifragility finding.

## Files

### Code

| File | Purpose |
|---|---|
| `run_verification.ps1` | Native C# implementation (compiled by Windows PowerShell 5.1) that checks the model's properties on randomized and edge cases and reproduces every reported value: worked example, one-at-a-time sensitivity, 2,000-draw global sensitivity, 30-period horizon, convexity/Jensen results, selection, feedback/bistability, and the capped/uncapped comparison. |
| `make_figures.py` | Python implementation of the worked-example model (mirrors `ExampleInp`, `Period`, `RunDesign` in the native script), Figures 1-4, the graphical abstract, `Cap_Comparison.csv` and `Illustrative_Trajectories.csv`. |
| `reproduce_scenario_analyses.py` | Recomputes all 2,000 native global-sensitivity draws in Python and produces the structural, matched-benchmark, dossier and scenario outputs. |
| `reproduce_numerical_checks.py` | Collects the native safeguard and boundary checks, replays the feedback maps, and produces the sampling-design, projection and summary outputs. |
| `test_model_outputs.py`, `test_numerical_checks.py`, `test_boundary_checks.py` | `pytest`/`unittest` suites (30 + 39 + 20 tests). |
| `requirements.txt` | Pinned dependencies (tested with Python 3.12.10). |

### Outputs

| File | Produced by | Content |
|---|---|---|
| `Verified_Native_Output.txt` | `run_verification.ps1` | Full native log: environment, script SHA-256, seeds, test denominators, all reported values, the 2,000 global-sensitivity draws (`GD|` rows) and feedback exports (`FEEDBACK|` rows). |
| `Cap_Comparison.csv` | `make_figures.py` | Five-year review with the saturation rule `min(1.5, annual_signal * review_years)` vs uncapped aggregation (disturbance 1.5 vs 4.0); both breach the group-harm limit. |
| `Illustrative_Trajectories.csv` | `make_figures.py` | Two governance designs x 31 timepoints; all six actual/counterfactual components, floors, performance and horizon flags. |
| `Global_Sensitivity_Draws.csv` | `reproduce_scenario_analyses.py` | The 2,000 native draws (25 columns) checked against Python. |
| `Model_Cross_Implementation_Checks.json` | `reproduce_scenario_analyses.py` | Python/C# agreement (8,000 numeric, 4,000 class-sequence, 8,000 horizon comparisons) and draw counts, e.g. capability gain at all five reviews in 358/2,000 draws (17.9%). |
| `Structural_Robustness.csv` | `reproduce_scenario_analyses.py` | 13 structural variants: routine-learning contrasts (0.120 at zero routine input to 0.038 at 0.8), monitoring/salience forms (0.103 baseline; geometric 0.169, minimum 0.129, capped-linear 0.097), coupled resources and reversibility floors. |
| `Matched_Decision_Benchmarks.csv` | `reproduce_scenario_analyses.py` | Five constructed benchmarks assessed with the same records and resources. |
| `Illustrative_Dossier.csv` | `reproduce_scenario_analyses.py` | Fictional programme: integrated burden 0.24 per year against a cumulative budget of 1.00; rejected in year 5 despite a positive capability contrast. |
| `Scenario_Checks.json` | `run_verification.ps1 -ScenarioChecksOnly` | 31 scenario checks (e.g. group harm 1.5 against a group limit of 1 with aggregate 0.75, temporary and cumulative protections, critical floors, routine learning, canonical evidence); nine are also checked against exact rational values in Python. |
| `Safeguard_Checks.json` | `run_verification.ps1 -SafeguardChecksOnly` | 87 checks of evidence fusion, numerical stability, risk tolerances, fixed-baseline horizon protection, required reporting routes and invalid-input rejection. |
| `Boundary_Checks.json` | `run_verification.ps1 -BoundaryChecksOnly` | 21 finite-boundary, rate-conversion, certain-harm, Hill-precision and qualified-theory checks. |
| `Feedback_Sampled_Maps.json`, `Feedback_Jacobians_and_States.json`, `Feedback_Local_Checks.csv` | `reproduce_numerical_checks.py` | Five sampled feedback systems, 25 replayed initial/terminal states, Jacobians/eigenvalues and 300 local perturbation returns (capacity subsystem only; not global basins). |
| `Sampling_Design_Stress.csv`, `Sampling_Design_Draws.csv` | `reproduce_numerical_checks.py` | Four artificial 512-point dependence/resource designs (not empirical priors). |
| `Projection_Diagnostics.csv` | `reproduce_numerical_checks.py` | Raw vs projected states; clipping an endpoint is not observed path safety. |
| `Numerical_Checks_Summary.json` | `reproduce_numerical_checks.py` | Summary: feedback checks, analytic bistability example, eta-q non-identifiability, cross-language boundary checks and sampling-design counts. |

## Running

Windows PowerShell, from this folder:

```powershell
python -m pip install -r requirements.txt pytest
python -m pytest -q -p no:cacheprovider                         # 89 passed
powershell -NoProfile -ExecutionPolicy Bypass -File .\run_verification.ps1   # RESULT: PASS, exit code 0
```

Optional switches of `run_verification.ps1` run a subset only:
`-InterfaceChecksOnly` (evidence interface, cap comparison and domain checks),
`-ScenarioChecksOnly`, `-SafeguardChecksOnly`, `-BoundaryChecksOnly` (each
prints its checks as JSON), and `-DefinitionsOnly` (compile only).

Regenerate the outputs into separate folders and compare them with the shipped
files:

```powershell
python -B make_figures.py .\generated
python -B reproduce_scenario_analyses.py --output .\recomputed
python -B reproduce_numerical_checks.py --output .\recomputed
```

The two reproduction drivers read `Verified_Native_Output.txt` and require the
SHA-256 recorded in it to match `run_verification.ps1`. If you modify the
script, rerun the full verification and save its output first, e.g.
`powershell -NoProfile -ExecutionPolicy Bypass -File .\run_verification.ps1 | Out-File -Encoding utf8 .\Verified_Native_Output.txt`.

On macOS/Linux, `make_figures.py` and the tests that do not call PowerShell
can be run (see the top-level README); the native script was only verified
under Windows PowerShell 5.1.

## Notes for interpretation

- The operational class is `G`/`CapabilityGain` under the Protected Attribution
  Rule, not volatility antifragility. Both arms receive the same valid routine
  input, so routine learning earns no disturbance credit.
- Index retention allows component trade-offs above predeclared critical
  floors. The reversibility floor of 0.30 is illustrative, not a calibrated safe
  level. Floors reject outcomes; they do not clip stocks upward.
- A protected horizon requires initial performance and capacity floors and, at
  every review, performance no lower than the fixed horizon-start baseline minus
  the declared tolerance (0.01 here). Required group reporting routes have an
  explicit veto under every monitoring aggregation.
- The fixed-direction Jensen/convexity result (inflection 0.231 for h = 2) needs
  a fixed state, zero routine input, valid noise-free evidence, bounded support
  below the inflection and acquisition dominating consumption. Radial convexity
  does not imply convexity in compound directions; the native script includes a
  counterexample.
- Feedback and bistability results are local numerical checks of the capacity
  subsystem, not full-state attraction or global basins.
- Figures 1-3 are conceptual diagrams; Figure 4 shows synthetic trajectories,
  not measurements from the five cases.
- The field `violations` on Observation/Numerical log rows is a reported count,
  not a failed assertion; test failures apply to Property, Implementation,
  Consistency and Edge rows.
