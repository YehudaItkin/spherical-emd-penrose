# Proof of the Einstein–Maxwell–dilaton Penrose inequality in spherical symmetry at α = 1

Source, manuscript, and machine-verifiable supplementary material for the paper

> **Proof of the Einstein–Maxwell–dilaton Penrose inequality in spherical symmetry at α = 1**
> Igor Itkin (independent researcher, [ig.itkin@gmail.com](mailto:ig.itkin@gmail.com))

## Summary

For the Einstein–Maxwell–dilaton (EMD) system at dilaton coupling α = 1, the paper proves the sharp
Pythagorean Penrose bound

$$2M \ge \sqrt{\rho_H^2 + 2q^2},$$

conjectured by Khuri, Weinstein and Yamada, for static, spherically symmetric, time-symmetric data
satisfying the dominant energy condition, with equality if and only if the data is the
Gibbons–Maeda / Garfinkle–Horowitz–Strominger (GHS) black hole. The paper also shows that the
divergence-free (source-free) hypothesis on the Maxwell field is *necessary*, while the
sub-extremality condition ρ_H ≥ q usually attached to charged Penrose inequalities is *redundant* here.

The bound is proved twice, by two monotone quasilocal masses of different algebraic type — a radical
Hawking–Bray mass (Section 5) and a radical-free rational mass (Section 7) — so each proof checks the
other. The second proof and Appendix A reduce to explicit polynomial positivity, which the
supplementary material verifies in exact arithmetic.

## Contents

| Path | Description |
|------|-------------|
| `spherical_emd_penrose.tex` | LaTeX source |
| `spherical_emd_penrose.pdf` | Compiled manuscript |
| `supplementary/` | Machine-verifiable certificate material for Section 7 and Appendix A |
| `supplementary/verify.py` | Self-contained sympy script re-deriving every algebraic claim in exact arithmetic |
| `supplementary/B9.txt`, `M9.txt` | The two polynomials whose nonnegativity is the monotonicity condition |
| `supplementary/M9_handelman_coefficients.txt` | The identity (1+t)⁸·M₉ written out: 1844 nonnegative coefficients |
| `supplementary/README.txt` | Details of the supplementary bundle |

## Build the manuscript

```bash
latexmk -pdf spherical_emd_penrose.tex
```

(No bibtex/biber step: the bibliography is an inline `thebibliography` environment.)

## Verify the computer-assisted steps

The proofs are analytic except for two polynomial-positivity certificates, which are exact over the
rationals (no floating point). To re-check them:

```bash
cd supplementary
python3 verify.py     # requires Python 3 with sympy
```

Expected final line: `ALL CHECKS PASSED.` The script verifies the Bernstein decomposition of B₉, the
one-line sum-of-squares for its top coefficient, the discriminant identities, the strict positivity of
the cofactors p₁₆, p₂₀ on (0, ∞) by Sturm's theorem, and the all-nonnegative Handelman certificate for M₉.

## Citation

A preprint is posted on Research Square. Please cite the manuscript; a BibTeX entry will be added here
once a DOI is assigned.

## License

This repository is dual-licensed:

- **Manuscript** — the paper text and figures (`spherical_emd_penrose.tex`, `spherical_emd_penrose.pdf`) are
  licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-CC-BY-4.0.txt).
- **Code** — the supplementary scripts (`supplementary/`, including `verify.py`) are licensed under the
  [MIT License](LICENSE).

© 2026 Igor Itkin.
