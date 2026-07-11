Supplementary material for

    "The spherically symmetric Einstein-Maxwell-dilaton Penrose inequality
     at alpha=1, and the redundancy of the charge-horizon hypothesis"
    Igor Itkin  (ig.itkin@gmail.com)

This bundle supports the second, radical-free proof (Section 7) and Appendix A.
Every object is exact over the rationals; nothing depends on floating point.

Variables:  t = e^phi > 0,   u = 2m/rho in (0,1),   s = sigma > 1.
Domain:     D = { t > 0, 0 < u < 1, s > 1 }.

Contents
--------
B9.txt
    The polynomial B9(t,u,s) whose nonnegativity on D is one of the two
    monotonicity conditions (num(beta_9) = t * B9). Text, ^ for powers.

M9.txt
    The polynomial M9(t,u,s) (num(Delta_9) = (u-1) t P^2 M9); its
    nonnegativity on D is the other monotonicity condition.

M9_handelman_coefficients.txt
    The identity eq:handelman of Section 7.3 written out in full:
        (1+t)^8 * M9 = sum_{k=0}^{3} Bcoef[k](t,w) * B_{3,k}(u),   w = sigma-1,
    with B_{3,k}(u) = C(3,k) u^k (1-u)^{3-k} >= 0 on [0,1]. All 1844
    coefficients of the four Bcoef[k](t,w) are listed and are nonnegative,
    so M9 >= 0 on D by inspection. This is the object the paper refers to as
    "provided as an ancillary file".

verify.py
    A self-contained sympy script that re-derives, in exact arithmetic, every
    algebraic claim of Section 7.3 and Appendix A from B9.txt and M9.txt:
      - the Bernstein-in-u decomposition of B9;
      - b2 = (t^2+1)^4 (s t+1)^2 r2 and the one-line SOS eq:r2sos;
      - the footnote identity r2|_{s=1} = t^2(t+1)^2((tau-3)^2+7), tau=t+1/t;
      - the discriminant identities eq:discid, and c0,c2,d0,d2 > 0;
      - p16, p20 strictly positive on (0,oo) (Appendix A, Sturm);
      - the all-nonnegative Handelman certificate for (1+t)^8 M9.
    It prints PASS/FAIL for each and exits nonzero on any failure.

Requirements and use
--------------------
    Python 3 with sympy.   Run:   python3 verify.py
    Expected final line:  ALL CHECKS PASSED.

The first proof (Section 5, the radical Hawking-Bray mass) is certified
separately by the exact Bernstein/Sturm argument described there and cross-
checked by cylindrical algebraic decomposition; it is independent of this bundle.
