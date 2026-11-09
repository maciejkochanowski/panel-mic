# Installation and verification

The library contains reusable deterministic calculations. The separate scientific package contains simulation counts, immutable generation provenance, corrected analysis provenance, data adapters, plots and article sources. Installing the library does not download research data or rerun simulations.

## Environments

The declared Python range begins at 3.11. CI is configured for Python 3.11, 3.12 and 3.13 on Windows and Linux. Check the release's actual CI results before treating a platform as verified. The minimum-dependency job uses Python 3.11 with NumPy 1.26.0 and SciPy 1.11.0; newer Python versions require compatible wheel versions.

The article environment uses Python 3.12.14, NumPy 2.3.5 and SciPy 1.18.1. requirements-article.txt pins the core numerical dependencies; the scientific archive retains its full operating-system-specific lock. Minimum-dependency and article-environment tests answer different compatibility questions.

    python -m pip install -r docs/library/requirements-article.txt
    python -m pip install .
    python examples/library/tutorial.py

## Tests and builds

In the lightweight source archive:

    python -m pip install ".[dev]"
    python -m pytest
    python -m build

The scientific checkout also includes workflow-dependent tests; use its reproduction commands for the complete scientific checks. Release packaging copies library tests explicitly and excludes tests depending on omitted scientific workflows or datasets.

Install each wheel or source archive into a new environment without the development source on PYTHONPATH. Run the CLI version command, tutorial, QC example and operation-schema checks there. Record artifact hashes and test receipts. A source-checkout test is not a test of a published wheel.

## Numerical and scientific boundaries

Conservative endpoints can differ from earlier versions by a few floating-point units. Decisions near a numerical boundary follow the certified direction, not an arbitrary epsilon. Unresolved numerical cases remain conservative or explicitly unavailable.

The simulation generation campaign is immutable. Reanalysis with corrected software preserves generation receipts and identifies the new analysis separately. Tests, hashes and replay do not establish IID sampling, biological transport, prospective validity or author approval.

## Publishing

The supplied publishing workflow is manually triggered. It builds and verifies distributions, then uses PyPI Trusted Publishing through a configured GitHub environment. Configure the repository, environment and trusted-publisher relationship in TestPyPI or PyPI before running it. No credentials are stored in source.

The workflow does not publish to Zenodo or assign a DOI. Update citation metadata only with an observed release URL and DOI for the actual artifact. Preserve exact versions and downloadable archive hashes.
