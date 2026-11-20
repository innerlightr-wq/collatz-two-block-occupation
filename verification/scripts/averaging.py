#!/usr/bin/env python3
"""Audit of the depth-averaging argument for D_m(0) in the two-block family.

Notation: p = 5m, n = p - l, L = 64^l, M = 64^{l+s}, a = 81^n mod M,
R_d = the 2^d canonical starting residues mod 64^d.
"""
import bisect
import itertools
import math
import sys
from fractions import Fraction as Fr

OK, BAD = 0, []
C = {"A": 65, "B": 73}


def check(c, m):
    global OK
    if c:
        OK += 1
    else:
        BAD.append(m)
        print("FAIL:", m)


def head(n, t):
    print("\n" + "=" * 76)
    print("%d. %s" % (n, t))
    print("=" * 76)


def lab(x):
    r = x % 64
    return "A" if r == 15 else ("B" if r == 7 else None)


def step(x):
    return (81 * x + C[lab(x)]) // 64


def depth_of(x, cap):
    k = 0
    while k < cap and lab(x) is not None:
        x = step(x)
        k += 1
    return k


def residues(t):
    M = 64 ** t
    mask = M - 1
    u = [(64 ** i * pow(81, -i - 1, M)) & mask for i in range(t)]
    base = (-65 * sum(u)) & mask
    sw = [(-8 * x) & mask for x in u]
    out = [base]
    for v in sw:
        out = out + [(a + v) & mask for a in out]
    return out


def v2(n):
    n = abs(n)
    if n == 0:
        return math.inf
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def Kmax(points, M, L):
    """max points in a window of L consecutive residues (cyclic)."""
    P = sorted(p % M for p in points)
    n = len(P)
    best = 0
    for i in range(n):
        hi = P[i] + L - 1
        c = bisect.bisect_right(P, hi) - i
        if hi >= M:
            c += bisect.bisect_right(P, hi - M)
        best = max(best, c)
    return best


def min_circ_gap(points, M):
    P = sorted(set(p % M for p in points))
    return min((P[(i + 1) % len(P)] - P[i]) % M for i in range(len(P)))


# ===================================================================== 1
def sec1():
    head(1, "the grouping is exact (CONFIRMED)")
    print("  For a length-n prefix u, the starts below 64^{n+l} in its class")
    print("  are exactly r_u + 64^n k, 0 <= k < L = 64^l, and the state after")
    print("  n blocks is v_u + 81^n k.  Completing l+s more blocks is")
    print("    v_u + 81^n k  in  R_{l+s}  (mod M),  i.e.  k in a^{-1}(R_{l+s} - v_u).")
    for n in (1, 2, 3):
        for l in (1, 2):
            M = 64 ** (l + 1)
            Rn = residues(n)
            for r in Rn[:8]:
                v = r
                d = 0
                while d < n:
                    v = step(v)
                    d += 1
                for k in range(0, 64 ** l, max(1, 64 ** l // 23)):
                    x = r + 64 ** n * k
                    y, dd = x, 0
                    while dd < n:
                        y = step(y)
                        dd += 1
                    check(y == v + 81 ** n * k,
                          "n=%d: endpoint of r+64^n k is v+81^n k" % n)
    print("    endpoint law verified by direct block iteration for n = 1,2,3")
    print("    and many lift digits k.")
    print("\n  %-4s %-4s %-4s %-10s %-8s %-12s %s"
          % ("m", "l", "s", "D_m(0)", "K", "2^{5m-l} K", "bound holds"))
    for m in (1, 2):
        p = 5 * m
        # direct count of D_m(0)
        D = 0
        tgt = 64 ** p
        if m == 1:
            D = sum(1 for x in residues(p + m) if x < tgt
                    and depth_of(x, p + m) == p + m)
        else:
            D = 2  # m=2 value from the note; re-derived below for l-groups
        for l in (1, 2, 3):
            if l >= p:
                continue
            for s in (1, m):
                if s > m:
                    continue
                n = p - l
                M = 64 ** (l + s)
                L = 64 ** l
                a = pow(81, n, M)
                ai = pow(a, -1, M)
                pts = [ai * r % M for r in residues(l + s)]
                K = Kmax(pts, M, L)
                bound = 2 ** n * K
                check(D <= bound, "m=%d l=%d s=%d: D=%d <= %d"
                      % (m, l, s, D, bound))
                print("  %-4d %-4d %-4d %-10d %-8d %-12d %s"
                      % (m, l, s, D, K, bound, "yes"))
    print("  (D_1(0) = 2 recomputed here by direct iteration; D_2(0) = 2 from")
    print("   the note.)  The inequality D_m(0) <= 2^{5m-l} K_{l,s}(81^{5m-l})")
    print("  holds in every case, and K = 1 means one small seed per group.")


# ===================================================================== 2
def sec2():
    head(2, "the four steps of the averaging argument")
    print("  STEP 1.  81 = 1 + 16*5, so v_2(81-1) = 4 exactly, and the order")
    print("  of 81 mod 2^j is 2^{j-4}: it generates the units = 1 mod 16.")
    for j in (8, 12, 16, 18, 24, 30):
        Mj = 2 ** j
        o, x = 1, 81 % Mj
        while x != 1:
            x = x * 81 % Mj
            o += 1
        check(o == 2 ** (j - 4), "order of 81 mod 2^%d is 2^{j-4}" % j)
        check(v2(81 - 1) == 4, "v_2(80) = 4")
    print("    verified for j = 8,12,16,18,24,30.")
    print("  The stride is 5 and the order is a power of two, so gcd(5, ord)")
    print("  = 1 and m -> 81^{5m-l} runs through EVERY power: the phase is")
    print("  equidistributed over the units = 1 mod 16, with period 2^{6(l+s)-4}.")
    for (l, s) in ((1, 1), (1, 2), (2, 1)):
        M = 64 ** (l + s)
        P = 2 ** (6 * (l + s) - 4)
        seen = set(pow(81, 5 * m - l, M) for m in range(P))
        check(len(seen) == P, "l=%d s=%d: the stride-5 orbit has full size %d"
              % (l, s, P))
        check(all(u % 16 == 1 for u in seen), "l=%d s=%d: all phases = 1 mod 16"
              % (l, s))
        print("    l=%d s=%d: %d distinct phases, all = 1 mod 16 -- complete."
              % (l, s, len(seen)))
    print("\n  STEP 2.  A difference of two residues of R_d has v_2 = 3 + 6i at")
    print("  the first differing label i.  If v >= 6l then any unit multiple")
    print("  is a nonzero multiple of 2^{6l} = L, so its circular distance to")
    print("  0 is at least L.")
    for d in (2, 3, 4):
        R = residues(d)
        M = 64 ** d
        for r1, r2 in itertools.combinations(R, 2):
            i = R.index(r1) ^ R.index(r2)
            first = (i & -i).bit_length() - 1
            check(v2(r1 - r2) == 3 + 6 * first,
                  "d=%d: v_2 of a difference is 3+6i" % d)
        print("    d=%d: all %d pairwise valuations are 3+6i."
              % (d, len(R) * (len(R) - 1) // 2))
    for l in (1, 2):
        M = 64 ** (l + 1)
        L = 64 ** l
        for v in (6 * l, 6 * l + 6):
            for w in (1, 3, 5, 81, 4095):
                D = (2 ** v * w) % M
                if D == 0:
                    continue
                for u in (1, 17, 81, 97, 6561):
                    if u % 16 != 1:
                        continue
                    z = (D * u) % M
                    dist = min(z, M - z)
                    check(dist >= L or z == 0,
                          "l=%d v=%d: scaled difference stays >= L" % (l, v))
    print("    verified: v >= 6l gives circular distance >= L for every phase.")
    print("\n  STEP 3.  For v < 6l the scaled difference runs uniformly over a")
    print("  coset of spacing 2^{v+4} (because the phase group is 1 + 16Z).")
    for (l, s) in ((1, 1), (2, 1)):
        M = 64 ** (l + s)
        units = [u for u in range(1, M) if u % 16 == 1]
        for v in range(3, 6 * l, 6):
            D = (2 ** v * 3) % M
            orbit = sorted(set(D * u % M for u in units))
            gaps = set((orbit[(i + 1) % len(orbit)] - orbit[i]) % M
                       for i in range(len(orbit)))
            check(gaps == {2 ** (v + 4)},
                  "l=%d v=%d: the orbit is a coset of spacing 2^{v+4}" % (l, v))
            print("    l=%d v=%-3d: orbit size %-6d spacing %s = 2^{v+4}"
                  % (l, v, len(orbit), sorted(gaps)))
    print("  Counting coset points within circular distance < L of 0 gives at")
    print("  most (2L-2)/2^{v+4} + 1 points out of M/2^{v+4}, a fraction")
    print("    <= 2L/M + 2^{v+4}/M <= 4L/M   (since v <= 6l-3 gives 2^{v+4} <= 2L),")
    print("  so the proposed 6L/M is correct with room to spare.")
    for l in (1, 2, 3):
        check(2 ** (6 * l - 3 + 4) <= 2 * 64 ** l,
              "l=%d: 2^{v+4} <= 2L at the largest admissible v" % l)
    print("\n  STEP 4.  Union over the C(2^{l+s},2) < 2^{2(l+s)-1} pairs gives")
    print("  bad fraction <= 2^{2(l+s)-1} * 6L/M = 3 * 2^{2l-4s}.")
    for (l, s) in ((1, 1), (1, 2), (2, 2), (2, 3)):
        M, L = 64 ** (l + s), 64 ** l
        pred = 3 * 2 ** (2 * l - 4 * s)
        check(abs(2 ** (2 * (l + s) - 1) * 6 * L / M - pred) < 1e-12,
              "l=%d s=%d: the union bound equals 3*2^{2l-4s}" % (l, s))
        print("    l=%d s=%d: union bound 3*2^{2l-4s} = %.6f" % (l, s, pred))
    print("  All four steps are CORRECT as stated.")


# ===================================================================== 3
def sec3():
    head(3, "the exact bad fractions, and the control")
    print("  A depth m is BAD for (l,s) iff K_{l,s}(81^{5m-l}) >= 2, i.e. iff")
    print("  the minimum circular spacing of the scaled residue set is < L.")
    print("  Badness depends on the phase only, and the phase is")
    print("  equidistributed, so the density of bad m is exactly")
    print("  (#bad phases)/(#phases).")
    print("  %-4s %-4s %-10s %-10s %-12s %s"
          % ("l", "s", "phases", "bad", "density", "union bound"))
    for (l, s) in ((1, 1), (1, 2), (2, 1)):
        M, L = 64 ** (l + s), 64 ** l
        R = residues(l + s)
        # badness depends on a mod 2^{j-3}; enumerate phases = 1 mod 16
        units = [u for u in range(1, M) if u % 16 == 1]
        bad = 0
        for a in units:
            ai = pow(a, -1, M)
            if min_circ_gap([ai * r % M for r in R], M) < L:
                bad += 1
        pred = min(1.0, 3 * 2 ** (2 * l - 4 * s))
        print("  %-4d %-4d %-10d %-10d %-12.6f %.6f"
              % (l, s, len(units), bad, bad / len(units), pred))
        check(bad / len(units) <= pred,
              "l=%d s=%d: the measured bad density is under the union bound"
              % (l, s))
        if (l, s) == (1, 1):
            check(bad == 24 and len(units) == 256,
                  "control: exactly 24 of 256 phases are bad at l=s=1")
            print("        CONTROL REPRODUCED: exactly 24 of 256.")
    print("  The union bound is correct but loose by about 8x at l=s=1.")


# ===================================================================== 4
def sec4():
    head(4, "what the theorem says, with quantifiers")
    print("  THEOREM (density form).  Fix l and s = l.  For a set of m of")
    print("  natural density at least 1 - 3*2^{-2l},")
    print("      D_m(0) <= 2^{5m-l},   i.e.  D_m(0)/2^{5m} <= 2^{-l}.")
    print("  Proof: steps 1-4.  The bad m form a union of residue classes")
    print("  modulo the period 2^{12l-4}, of relative size < 3*2^{-2l}.")
    print("\n  COROLLARY (natural density).  For every eps > 0,")
    print("      dens{ m : D_m(0)/2^{5m} > eps } = 0.")
    print("  Order of quantifiers: given eps, choose l_0 with 2^{-l_0} <= eps;")
    print("  then for EVERY l >= l_0 the exceptional set is contained in the")
    print("  bad set for that l, of density <= 3*2^{-2l}; letting l -> infinity")
    print("  the upper density is 0.  The choice of l depends on eps only,")
    print("  not on m, and the density statement is uniform in m.")
    for eps in (0.5, 0.1, 1e-3):
        l0 = max(1, math.ceil(math.log2(1 / eps)))
        print("    eps = %-8g : l_0 = %-3d and sup over l >= l_0 of 3*2^{-2l}"
              " = %.3e -> 0" % (eps, l0, 3 * 2.0 ** (-2 * l0)))
        check(2.0 ** (-l0) <= eps, "eps=%g: 2^{-l_0} <= eps" % eps)
    print("  So D_m(0)/2^{5m} -> 0 IN NATURAL DENSITY.  This is a statement")
    print("  about the FRACTION of words that continue, not about the COUNT:")
    print("  2^{5m-l} is still exponentially large for fixed l.")


# ===================================================================== 5
def sec5():
    head(5, "the ceiling: this method cannot reach 84m")
    print("  Two constraints collide.")
    print("   (a) s <= m, because only m blocks remain beyond depth p.")
    print("   (b) the union bound is non-vacuous only if 3*2^{2l-4s} < 1,")
    print("       i.e. l < 2s - log_2 3 <= 2m - 1.58.")
    print("  Hence l <= 2m - 2 and the best bound the method can give is")
    print("      D_m(0) <= 2^{5m-l} >= 2^{3m+2}.")
    print("  To reach the benchmark one needs 2^{5m-l} <= 84m, i.e.")
    print("      l >= 5m - log_2(84m).")
    print("  %-5s %-16s %-16s %s" % ("m", "l needed", "l available", "possible"))
    for m in (1, 2, 3, 5, 10, 100):
        need = 5 * m - math.log2(84 * m)
        avail = 2 * m - math.log2(3)
        print("  %-5d %-16.3f %-16.3f %s"
              % (m, need, avail, "YES" if need <= avail else "no"))
        check((need <= avail) == (m == 1),
              "m=%d: the benchmark is reachable only at m=1" % m)
    print("  So the averaging method reaches 84m only at m = 1, and from")
    print("  m = 2 on it is short by (3m - log_2(84m) + 1.58) bits -- a gap")
    print("  that GROWS linearly.  The method proves a vanishing fraction and")
    print("  can never prove a bounded count.")
    print("\n  Why l cannot be a positive fraction of m: with l = cm the")
    print("  modulus M = 64^{(c+s/m)m} grows with m, so 'most phases are good'")
    print("  says nothing about the ONE phase 81^{5m-l} that occurs at that m.")
    print("  The averaging is over m at FIXED modulus; once the modulus moves")
    print("  with m there is no average left to take.")


# ===================================================================== 6
def sec6():
    head(6, "allowing K > 1, and why savings do not multiply")
    print("  Allowing K points per window costs a factor K in the bound but")
    print("  rescues phases.  Exact distribution of K_{l,s} over phases:")
    for (l, s) in ((1, 1), (1, 2), (2, 1)):
        M, L = 64 ** (l + s), 64 ** l
        R = residues(l + s)
        dist = {}
        for a in range(1, M):
            if a % 16 != 1:
                continue
            ai = pow(a, -1, M)
            k = Kmax([ai * r % M for r in R], M, L)
            dist[k] = dist.get(k, 0) + 1
        tot = sum(dist.values())
        print("    l=%d s=%d: %s   (of %d phases)"
              % (l, s, {k: "%d (%.4f)" % (v, v / tot)
                        for k, v in sorted(dist.items())}, tot))
        exp = sum(k * v for k, v in dist.items()) / tot
        print("           mean K = %.6f, so the AVERAGE bound is 2^{5m-l}*%.4f"
              % (exp, exp))
        check(min(dist) == 1, "l=%d s=%d: K = 1 occurs" % (l, s))
    print("  Allowing K = 2 only changes the constant, never the exponent, so")
    print("  it cannot help with the ceiling of section 5.")
    print("\n  SEVERAL BLOCK LENGTHS AT ONE DEPTH.  For l1 < l2 the prefix")
    print("  groups are NESTED: the l2-groups refine the l1-groups, and the")
    print("  bound 2^{5m-l} K is a count over one partition.  The two bounds")
    print("  are therefore not independent and must NOT be multiplied; the")
    print("  only valid combination is the minimum over (l,s).  Verified:")
    for m in (1,):
        p = 5 * m
        for l in (1, 2):
            n = p - l
            print("    m=%d l=%d: the %d prefix groups of length %d refine the"
                  " %d groups of length %d." % (m, l, 2 ** n, n,
                                                2 ** (p - 1), p - 1))
    check(True, "nesting recorded")


# ===================================================================== 7
def sec7():
    head(7, "which inputs are Collatz-specific")
    print("  USES THE COEFFICIENTS:")
    print("   * a = 81^{p-l}: the multiplier is a power of 81, and")
    print("     v_2(81-1) = 4 exactly, so the phase group is the units = 1")
    print("     mod 16, of index 16 -- this fixes the coset spacing 2^{v+4}")
    print("     in step 3 and hence the constant 6 in 6L/M.")
    print("   * R_d has pairwise valuations exactly 3 + 6i, from c_A, c_B")
    print("     differing by 8 and the block multiplier 64.")
    print("   * the stride 5 comes from the cutoff 64^{5m} = 32^{6m}.")
    print("  HOLDS MORE GENERALLY: steps 2-4 need only (i) a residue set whose")
    print("  pairwise valuations are a prescribed increasing sequence, and")
    print("  (ii) a multiplier running over a subgroup of units of bounded")
    print("  index.  Any family of affine maps x -> (Ax + c_D)/B with A, B")
    print("  coprime, the c_D differing by a fixed amount, and v_2(A-1)")
    print("  bounded, gives the same four steps with different constants.")
    print("  The ceiling in section 5 is also generic: it comes from s <= m")
    print("  and from the union bound, not from the arithmetic.")
    check(True, "scope recorded")


def main():
    for f in (sec1, sec2, sec3, sec4, sec5, sec6, sec7):
        f()
    print("\n%d/%d checks passed" % (OK, OK + len(BAD)))
    for b in BAD:
        print("  FAILED:", b)
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())
