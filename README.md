# Protected Attribution Rule: model and case records

This repository contains the code and data for a journal article on a
Protected Attribution Rule for sustainability assessment (citation details
will be added on publication). `model_verification/` holds the illustrative dynamic
learning-capacity model, the figure code, the native and Python verification
of the analytical and numerical results, and the outputs reported in the
article (Supplementary Material 1). `case_records/` holds the documentary
assessment of five Waterloo Region cases: integration ratings by two raters,
evidence-sufficiency coding on eight criteria, and the script that reproduces
the rating statistics and case statuses (Supplementary Material 2). File names
here are descriptive and may differ from those in the article's supplementary
archives; the code, data and results are the same. Please
cite the article once it is published; until then, see
[CITATION.cff](CITATION.cff).

## Scope

- The model is illustrative and uncalibrated. Synthetic simulation outputs are
  not empirical validation, and sampling frequencies are not real-world
  prevalence.
- The cases are documentary. No field measurements or interviews were
  conducted.
- Third-party source documents are cited with locators in the case records.
  They are not redistributed here.

## Folder map

```
.
|-- README.md              this file
|-- LICENSE                MIT licence (original code)
|-- LICENSE-DATA.md        CC BY 4.0 licence (original data and documentation)
|-- CITATION.cff           citation metadata
|-- .gitignore
|-- model_verification/    model, figures, verification scripts, tests and outputs (see its README)
`-- case_records/          case-record scripts and tests (see its README)
    |-- data/              coded inputs
    `-- results/           reproduced outputs
```

## Requirements

- Python 3.12 (tested with 3.12.10)
- `model_verification/requirements.txt`: `numpy==2.5.4`, `matplotlib==3.11.2`, `Pillow==12.3.0`
- `case_records/requirements.txt`: `numpy==2.5.4`
- `pytest` to run the test suites
- Windows PowerShell 5.1 (.NET Framework, included with Windows) for the native
  verification script and for the Python checks that call it. No .NET SDK is
  required. The figures use the Arial font.

## Reproduction

### Windows PowerShell

```powershell
cd model_verification
python -m pip install -r requirements.txt pytest
python -m pytest -q -p no:cacheprovider
powershell -NoProfile -ExecutionPolicy Bypass -File .\run_verification.ps1
cd ..

cd case_records
python -m pip install -r requirements.txt pytest
python -m pytest -q -p no:cacheprovider
python -B reproduce_cases.py --output .\recomputed
cd ..
```

### macOS / Linux

The native verification script and the tests that call it need Windows
PowerShell, so they are excluded here; the remaining tests and the case-record
reproduction run on any platform.

```bash
cd model_verification
python3 -m pip install -r requirements.txt pytest
python3 -m pytest -q -p no:cacheprovider --ignore=test_boundary_checks.py \
  --deselect test_model_outputs.py::ModelTests::test_native_scenario_checks_and_exact_vectors
cd ..

cd case_records
python3 -m pip install -r requirements.txt pytest
python3 -m pytest -q -p no:cacheprovider
python3 -B reproduce_cases.py --output ./recomputed
cd ..
```

Commands for regenerating every output file are listed in
[model_verification/README.md](model_verification/README.md) and
[case_records/README.md](case_records/README.md).

### Expected results

| Folder | Command | Expected |
|---|---|---|
| `model_verification` | `python -m pytest -q` (Windows) | 89 passed |
| `model_verification` | macOS/Linux command above | 68 passed, 1 deselected |
| `model_verification` | `run_verification.ps1` | `RESULT: PASS`, exit code 0 (about 1-2 minutes) |
| `case_records` | `python -m pytest -q` | 24 passed |

## Licences

- Original code: MIT - see [LICENSE](LICENSE).
- Original data and documentation: CC BY 4.0 - see [LICENSE-DATA.md](LICENSE-DATA.md).
- Third-party documents are not redistributed and remain under their owners' terms.

## Citation

Please cite the article once it is published. Until then, cite this code and
data using the metadata in [CITATION.cff](CITATION.cff).
