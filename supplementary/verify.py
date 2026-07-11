#!/usr/bin/env python3
"""
Supplementary verification for

    "The spherically symmetric Einstein-Maxwell-dilaton Penrose inequality
     at alpha=1, and the redundancy of the charge-horizon hypothesis"
    (Igor Itkin).

This script re-derives, in EXACT rational arithmetic (sympy), every algebraic
claim of Section 7.3 ("The rational certificate") and Appendix A, from the two
polynomials B9 and M9 supplied in B9.txt and M9.txt.  It prints PASS/FAIL for
each check and exits nonzero on any failure.

Variables:  t = e^phi > 0,  u = 2m/rho in (0,1),  s = sigma > 1
Domain:     D = { t>0, 0<u<1, s>1 }.

Requires: sympy.  Run:  python3 verify.py
"""
import sys
import sympy as sp

t, u, s, w, tau = sp.symbols('t u s w tau')


def load(fname):
    with open(fname) as f:
        return sp.expand(sp.sympify(f.read().replace('^', '**'),
                                    locals={'t': t, 'u': u, 's': s}))


def check(name, ok):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    return bool(ok)


def positive_on_open_halfline(P):
    """True iff the univariate polynomial P(t) is > 0 on (0, oo)."""
    P = sp.expand(P)
    poly = sp.Poly(P, t)
    m = min(mon[0] for mon in poly.monoms())          # strip t^m (root at 0)
    core = sp.Poly(sp.expand(P / t**m), t)
    return core.count_roots(0, sp.oo) == 0 and core.eval(sp.Integer(1)) > 0


def main():
    B9 = load('B9.txt')
    M9 = load('M9.txt')
    results = []

    print("B9 >= 0 on D  (Section 7.3, eq. for B9):")
    # Bernstein-in-u:  B9 = b0 (1-u)^2 + 2 b1 u(1-u) + b2 u^2
    c = sp.Poly(B9, u).all_coeffs()[::-1]              # a0 + a1 u + a2 u^2
    a0, a1, a2 = c[0], c[1], c[2]
    b0 = sp.expand(a0)
    b2 = sp.expand(a0 + a1 + a2)
    b1 = sp.expand(a0 + a1 / 2)
    recon = sp.expand(b0 * (1 - u)**2 + 2 * b1 * u * (1 - u) + b2 * u**2)
    results.append(check("Bernstein-in-u decomposition of B9", sp.expand(recon - B9) == 0))

    # b2 = (t^2+1)^4 (s t+1)^2 r2  with one-line SOS (eq:r2sos)
    r2 = sp.cancel(b2 / ((t**2 + 1)**4 * (s * t + 1)**2))
    results.append(check("b2 = (t^2+1)^4 (s t+1)^2 r2, r2 polynomial",
                         sp.denom(sp.cancel(r2)) == 1))
    C2 = t**2 * ((t**2 - 2)**2 + 7)
    C1 = -4 * t * (t**2 - 2 * t - 1) * (t**2 + 2 * t - 1)
    sos = sp.expand((2 * C2 * s + C1)**2 + 28 * t**2 * (t**2 - 1)**2 * (t**2 + 1)**2)
    results.append(check("4 C2 r2 = (2 C2 s + C1)^2 + 28 t^2(t^2-1)^2(t^2+1)^2 (eq:r2sos)",
                         sp.expand(4 * C2 * sp.expand(r2) - sos) == 0))
    results.append(check("C2 = t^2((t^2-2)^2+7) > 0 on t>0", positive_on_open_halfline(C2)))
    # r2|_{s=1} = t^2 (t+1)^2 ((tau-3)^2 + 7),  tau = t + 1/t   (Section 7.3 footnote)
    r2_s1 = sp.expand(r2.subs(s, 1))
    claim = t**2 * (t + 1)**2 * (((t + 1 / t) - 3)**2 + 7)
    results.append(check("r2|_{s=1} = t^2(t+1)^2((tau-3)^2+7)  (tau=t+1/t)",
                         sp.simplify(r2_s1 - claim) == 0))

    # b1 = (t^2+1)^2 (s t+1) r1 ;  discriminant identities (eq:discid)
    r1 = sp.expand(sp.cancel(b1 / ((t**2 + 1)**2 * (s * t + 1))))
    cc = sp.Poly(r1.subs(s, 1 + w), w).all_coeffs()[::-1]   # c0 + c1 w + c2 w^2 + ...
    dd = sp.Poly(b0.subs(s, 1 + w), w).all_coeffs()[::-1]   # d0 + d1 w + ...
    disc_c = sp.expand(cc[1]**2 - 4 * cc[0] * cc[2])
    disc_d = sp.expand(dd[1]**2 - 4 * dd[0] * dd[2])
    p16 = sp.expand(sp.cancel(-disc_c / (t**2 * (t + 1)**4)))
    p20 = sp.expand(sp.cancel(-disc_d / (4 * t**2 * (t + 1)**6 * (t**2 + 1)**2)))
    results.append(check("c1^2-4 c0 c2 = -t^2(t+1)^4 p16 (eq:discid)",
                         sp.expand(disc_c + t**2 * (t + 1)**4 * p16) == 0))
    results.append(check("d1^2-4 d0 d2 = -4 t^2(t+1)^6(t^2+1)^2 p20 (eq:discid)",
                         sp.expand(disc_d + 4 * t**2 * (t + 1)**6 * (t**2 + 1)**2 * p20) == 0))
    results.append(check("c0,c2 > 0 and higher (c3) > 0 on t>0",
                         all(positive_on_open_halfline(cc[i]) for i in (0, 2, 3))))
    results.append(check("d0,d2 > 0 and higher (d3,d4) > 0 on t>0",
                         all(positive_on_open_halfline(dd[i]) for i in (0, 2, 3, 4))))

    print("\nAppendix A  (p16, p20 strictly positive on (0,oo) by Sturm):")
    results.append(check("p16 has no root in (0,oo) and p16(1)=3072>0",
                         positive_on_open_halfline(p16) and p16.subs(t, 1) == 3072))
    results.append(check("p20 has no root in (0,oo) and p20(1)=8192>0",
                         positive_on_open_halfline(p20) and p20.subs(t, 1) == 8192))

    print("\nM9 >= 0 on D  (Section 7.3, Handelman identity eq:handelman):")
    # (1+t)^8 M9, with s = 1+w, expanded in the cubic u-Bernstein basis,
    # has all-nonnegative coefficients in (t,w).
    P = sp.expand(M9.subs(s, 1 + w) * (1 + t)**8)
    from sympy import binomial as C
    cs = sp.Poly(P, u).all_coeffs()[::-1]
    total = 0
    allnn = True
    for j in range(4):
        bj = sp.expand(sum(C(j, kk) * cs[kk] / C(3, kk) for kk in range(j + 1)))
        pj = sp.Poly(bj, t, w)
        total += len(pj.terms())
        if any(co < 0 for _, co in pj.terms()):
            allnn = False
    results.append(check(f"(1+t)^8 M9 in u-Bernstein basis: all {total} (t,w)-coefficients >= 0",
                         allnn))

    print()
    if all(results):
        print("ALL CHECKS PASSED.")
        return 0
    print("SOME CHECKS FAILED.")
    return 1


if __name__ == '__main__':
    sys.exit(main())
