# Exact finite L-ensemble contract

## Established method

A real symmetric positive-semidefinite matrix L indexed by n items defines an L-ensemble:

```
P(Y=S) = det(L[S]) / Z
Z = det(I+L) = sum over all S of det(L[S])
K = L(I+L)^-1 = I-(I+L)^-1
P(A subset Y) = det(K[A])
```

The empty determinant is one, so every finite L-ensemble assigns positive probability to empty output. This is a proper subset of all DPPs; projection inclusion kernels with eigenvalue one do not have finite L representations. We do not accept a K matrix as if it were L.

Native determinants use exact Fraction Gaussian elimination with row swaps; inverses use Gauss-Jordan. Every principal minor is checked nonnegative, a necessary and sufficient PSD criterion for real symmetric matrices. Leading minors alone are not enough for semidefiniteness. Singular PSD matrices are valid. There is no epsilon, clipping, jitter or eigenvalue repair. The complete bounded subset table is retained in size-then-lexicographic order, including zero-weight subsets.

## Inclusion, not exact-set probability

For distinct i,j, pair inclusion is Kii*Kjj-Kij*Kji and covariance is -Kij². The equal-marginal independent pair probability Kii*Kjj is retained as a control, not a competing fitted model. Either off-diagonal sign can imply repulsion. E|Y| = trace(K); Var|Y| = trace(K)-trace(K²). Enumeration independently checks every marginal and both moments.

## Include/exclude conditioning

Let A be required included, B required excluded, R the remaining indices. Drop B from L, then condition A using the Schur complement:

```
C = L[R,R] - L[R,A] inverse(L[A,A]) L[A,R]
```

For positive-probability events, L[A,A] is invertible, including the 0x0 inverse for empty A. C is the residual L-ensemble likelihood kernel on R; its inclusion kernel is I-(I+C)^-1. The event weight is the exact sum of original det(L[S]) for S including A and avoiding B. Conditional atom probabilities divide by that event weight. Direct enumeration verifies every residual atom and every inclusion minor. If event weight is zero, both residual matrices and the law are undefined/null. No pseudoinverse is used.

## Optional fixed total size

The k field, when set, additionally requires |Y|=k, including already-required items. Selected probabilities are direct determinant weights normalized over the intersection of all conditions. When no positive weight satisfies that event, selected laws/marginals/moments are null, modes empty. Zero probabilities for excluded atoms in a possible event are ordinary zero strings, not null.

The residual C and its K condition only on include/exclude. They do not incorporate k. Fixed-size k-DPPs are not ordinary DPPs in general; the final selected marginals and size law come from exact enumeration, not minors of the event-only K. Common positive scaling L -> cL multiplies each size-k determinant by c^k and cancels in the fixed-k law. It does not generally cancel in the variable-size law. No efficient k-DPP sampler or eigensolver is implemented.

## Modes and complete output

All exact maximum-weight subsets are reported, with ties retained. Modes use original weights; the normalizer is common. Selected modes maximize over the positive-probability selected event. A mode is not a draw. No greedy MAP approximation is implemented and no global efficient-MAP claim is made. No ranking or relevance quality was evaluated.

## Strict input and exactness boundary

The validated public entry point is audit(spec). Exactly L/include/exclude/k are required. Matrix size 1..5. Matrix cells must be strings matching an integer or integer/positive-integer fraction; no decimal strings, exponent notation, whitespace, Booleans or numeric JSON cells. Strings <=12 characters; reduced numerator magnitude and denominator <=10,000. Equivalent fractions normalize internally, while spec preserves supplied strings. Include/exclude indices are distinct integers in range, disjoint. k is null or an integer 0..n. Python Booleans and fractional values are refused for indices and k; finite integral numbers such as JSON 1.0 are accepted just as JS integral Number values are. Integer-valued indices normalize to integers before matrix addressing.

Native input is <=32KiB, rejects duplicate keys and nonfinite tokens. Browser input is <=32,768 characters; JSON.parse keeps the last duplicate key. This is not a byte-identical parser claim. All computation and report fractions are exact rational arithmetic. Browser SVG widths convert fractions to Number for geometry only; complete exact strings remain inspectable. No sampling, training, datasets or API call is involved. Runtime is exponential/factorial for verification and intentionally tiny.

Valid impossible events exit CLI 0 with null conditional results. Invalid input exits 2. No quality threshold is imposed.

## Independent checks

24 Python tests run on 3.10-3.13. Independent permutation expansion checks 120 determinants (including empty and nonsymmetric matrices); 50 inverse-product identities; 140 seeded Gram kernels check all subsets, marginal determinant identities, exact normalization and size moments; 100 independently normalized include/exclude/k controls check all selected atoms and event-only residual atom/inclusion laws; 40 common-scale fixed-k controls; 20 sign-conjugation controls. Diagonal independence, duplicate rank-one support, negative off-diagonal repulsion, impossible events, all-in/all-out and strict rejection remain explicit. These counts are control counts within tests, not separate test methods.

24 complete reports compare 4,050 exact leaves against independent BigInt JavaScript. JS uses permutation expansion and adjugate inverses, not the native elimination algorithm. Actual Chrome checks all 252 subset rows and full records, 786 matrix cells with full selected matrix records, 84 pair rows, 504 bar fraction/geometry records, selected marginal/cardinality/mode records, all cards, invalid/stale behavior, JSON export and 320/390/1280 containment. Isolated wheel reproduction runs outside the source checkout, checks all retained reports, installed console CLI and four installed module origins.

## Sources, read November 9, 2026

1. Kulesza and Taskar (2012/2013), [Determinantal point processes for machine learning](https://arxiv.org/html/1207.6083v4), primary survey, sections 2.1-2.4: L-ensemble normalization, inclusion kernel, negative correlations and include/exclude conditioning. High confidence; main mathematical source. No code or figures copied.
2. [DPPy definitions](https://dppy.readthedocs.io/en/stable/finite_dpps/definition.html), official project documentation: likelihood vs inclusion, sufficient kernel conditions and k-DPP distinction. High confidence in explicit formulas. Example-code comments/notation sometimes use marginal ambiguously; this project uses exact-set versus inclusion consistently.
3. [DPPy properties](https://dppy.readthedocs.io/en/latest/finite_dpps/properties.html), official project documentation: pair covariance, trace moments and conditional inclusion kernels. High confidence; no software dependency introduced.
4. Swati Gupta, MIT 6.883 Fall 2015, [Lecture 21](https://people.csail.mit.edu/stefje/fall15/notes_lecture21.pdf), independent educational notes: normalizer proof, kernel mapping, geometry and conditioning; corroborates primary formulas. PDF extraction repeats some slides and damages math glyphs; primary formulas control.
5. Gautier et al. (2019), [DPPy: DPP Sampling with Python](https://jmlr.org/papers/volume20/19-179/19-179.pdf), JMLR software paper: established exact/approximate DPP sampling ecosystem and scale limits. This project is an exact finite-law audit, not an implementation or performance comparison of DPPy's samplers.
6. Buntine (2017-09-14), [Auxiliary notes](https://bayesian-models.org/wp-content/uploads/2017/09/dppaux_140917.pdf), independent lecture slides: Gram/Schur interpretation and representation distinctions. Caveat: extracted volume wording and singleton exclusion signs are ambiguous or damaged. We use squared Gram volume and the primary L-restriction/Schur identities; no code, figures or questionable formula copied.
7. [Semidefiniteness criterion](https://mjo.osborne.economics.utoronto.ca/index.php/tutorial/index/1/qfs/t), University of Toronto mathematical methods notes, proposition 3.2.3.1: all principal minors, not only leading minors, characterize positive semidefiniteness for symmetric matrices.
8. [DPPy exact sampling](https://dppy.readthedocs.io/en/latest/finite_dpps/exact_sampling.html), official algorithm documentation for contrast: sampling is separate from enumerating a complete law. No sampler implemented or benchmarked.

Implementation, independent checks and original synthetic controls by Midya Rahmani. Established mathematics, not a new algorithm, trained model or downstream result.
