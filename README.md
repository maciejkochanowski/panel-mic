# Panel MIC

Panel-MIC is a Python library and command-line interface for exact, conditional certification of upper MIC quantiles and planning future certification with ambiguous observations. It accompanies the methodological article **Exact study planning for upper quantile certification with ambiguous MIC observations**.

The method distinguishes whether the current data support a specified P99 or P99.9 threshold from whether a future study can attain a specified certification power. Current certification can succeed while future planning remains unresolved. Unknown observations remain in the denominator.

## Installation

Use a virtual environment with Python 3.11–3.13:

    python -m venv .venv

Activate the environment, then install the supplied source directory:

    python -m pip install .
    panel-mic --version

Core requirements are NumPy ≥1.26 and SciPy ≥1.11; the resolver must select versions compatible with your Python interpreter. The article's pinned environment is separate from these minimum requirements. See [verification notes](docs/library/REPRODUCIBILITY.md). A configured CI matrix is not itself evidence that every platform has passed.

## Example

    from panel_mic import Reference, assess_upper, design_upper

    reference = Reference([1.0], [392, 0], iid=True, unit="mg/L")
    current = assess_upper(reference, 1.0, percentile=.99, alpha=.025)
    future = design_upper(
        reference, [1.0], future_n=874,
        percentile=.99, target_power=.8, alpha=.025, delta=.05,
    )
    print(current["status"])  # PANEL_REACHES_QUANTILE
    print(future["status"])   # REFERENCE_N_TOO_SMALL

These counts are a synthetic illustration. The IID declaration requires scientific justification; setting it to true does not establish independent or representative sampling.

Corresponding CLI examples:

    panel-mic assess-upper examples/library/current-certification.json --text
    panel-mic design-upper examples/library/future-design.json --text
    panel-mic diagnose-upper examples/library/future-diagnostic.json --text
    panel-mic protocol-assurance examples/library/protocol-assurance.json --text

JSON is the default output. The --text option selects English explanations; add --language pl for Polish. Use --output result.json to create a new file. Existing files are never overwritten. Input errors produce exit code 2; a valid but inconclusive scientific result has exit code 0. Batch output retains every row, including invalid inputs.

## Scope

The production interface supports P99 and P99.9 and the article's primary planning scope of future N ≤1000. Larger numerical searches are available where documented; they do not extend the biological evidence. Fixed-scenario selection error is an unconditional repeated-procedure bound, not a posterior probability that a particular recommendation has sufficient power.

The library implements deterministic certification, bounded population changes, future ambiguity budgets, conditional QC profiles, and fixed-protocol assurance. Randomization comparators and the full article reproduction workflow belong to the scientific supplement. No prospective biological validation or transport of QC profiles is implied.

## Documentation

- [English tutorial](docs/library/TUTORIAL.md): data, current certification, future power and interpretation.
- [Python API](docs/library/API.md): all exported operations and assumptions.
- [Input and output contract](docs/library/INPUT_CONTRACT.md): schemas, strict counts and status semantics.
- [Reproducibility](docs/library/REPRODUCIBILITY.md): installation, testing and environment boundaries.
- [Release changes](CHANGELOG.md) and [citation metadata](CITATION.cff).

The package uses the MIT license. Cite the exact software version and the methodological article when its final citation is available. Release metadata must identify the actual published archive; no DOI is implied by this source tree.
