# DPPLab

Exact finite determinantal point processes with every subset and conditional fraction exposed. Dependency-free Python, independent BigInt JavaScript, no model or API call.

[Live inspector](https://midya8.github.io/dpplab/) · [Method, sources and limits](METHOD.md)

## What it does

- Checks symmetric positive-semidefinite likelihood kernels by every exact principal minor. Valid singular kernels stay singular.
- Enumerates every subset with determinant weight and exact-set probability; compares every inclusion determinant against a separate atomic-probability sum.
- Computes K = I - (I+L)^-1, baseline pair covariances and cardinality moments.
- Conditions on included and excluded items, with explicit Schur-complement residual kernels and zero-event refusal.
- Adds optional fixed **total** size k by exact enumeration. Residual inclusion kernels deliberately remain event-only: k-DPPs are not ordinary DPPs in general.
- Retains every exact MAP tie, all selected item marginals, complete cardinality laws and JSON export.

The default control has two identical items. Each has inclusion probability 4/9, but their joint inclusion probability is zero, versus 16/81 from an independent equal-marginal control. A third independent item has inclusion probability 1/2. Supplied geometry, not measured semantic diversity.

## Run

```bash
python -m pip install .
dpplab --spec examples/duplicates.json --json report.json
# or: python -m dpplab --spec examples/duplicates.json --json report.json
```

```python
from dpplab import audit
report = audit({"L": [["2", "-1"], ["-1", "2"]],
                "include": [], "exclude": [], "k": None})
```

The exact root keys are L/include/exclude/k. One to five items, at most 32 subsets. Matrix entries are integer/fraction strings, at most 12 characters; reduced numerator magnitude and denominator are <=10,000. Indices are distinct zero-based integer-valued numbers, include and exclude disjoint. k is null or an integer 0..n. Native input is <=32KiB, rejects duplicate keys/nonfinite JSON tokens. Browser JSON.parse keeps the last duplicate value and bounds input characters, not UTF-8 bytes. Fractions normalize exactly; decimal approximation is display-only for bar widths.

A valid impossible event completes with null conditional laws/moments and exit 0. Invalid specifications exit 2. There is no quality gate, sampler or random seed API.

## Verification

```bash
python -m unittest discover -s tests -v
python tools/build_audits.py
npm test
bash tools/wheel_smoke.sh
```

24 Python tests pass on 3.10-3.13. They include 120 independent permutation determinants, 50 inverse identities, 140 Gram-kernel full subset/inclusion/size oracles, 100 conditioning and residual-kernel oracles, 40 fixed-size scaling controls and 20 sign-conjugation controls. Independent JavaScript uses permutation determinants and adjugate inverses, unlike native elimination, and checks 24 complete reports with 4,050 exact leaves. Actual Chrome checks 252 subset rows and full records, 786 matrix cells, 84 pair rows, 504 bar fractions/positions/widths, selection/mode records, all cards, stale/invalid input, export and 320/390/1280 containment. Isolated wheel reproduction checks all reports, console CLI and four installed module origins outside the checkout.

## Scope

Established real-symmetric finite L-ensembles only, not all DPPs. Full enumeration is exponential and intentionally tiny. No learned kernel, search corpus, relevance metric, sampler, large-scale inference, fairness result or downstream quality claim. Mode enumeration is not sampling. Exact rational results describe supplied synthetic mathematics, not floating-point embedding pipelines.

Implementation and synthetic controls by Midya Rahmani. No external implementation code, figures, datasets or model weights copied. Primary and secondary sources, numerical conventions and source contradictions are documented in [METHOD.md](METHOD.md).
