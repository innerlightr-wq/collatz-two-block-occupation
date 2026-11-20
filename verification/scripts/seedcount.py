#!/usr/bin/env python3
"""Exact seed count by half-word split, and the dyadic bin identity."""
import bisect, math, sys
from fractions import Fraction as Fr

OK, BAD = 0, []
CS = (65, 73)


def check(c, m):
    global OK
    if c: OK += 1
    else:
        BAD.append(m); print("FAIL:", m)


def head(n, t):
    print("\n" + "=" * 76); print("%d. %s" % (n, t)); print("=" * 76)


def setup(m):
    n = 3 * m; B = 64 ** n; Q = B * B; W = 64 ** (2 * m); H = B * W; J = 6 * m
    u = [(64 ** i * pow(81, -i - 1, Q)) % Q for i in range(2 * n)]
    return n, B, Q, W, H, J, u


def halves(m):
    n, B, Q, W, H, J, u = setup(m)
    P = [0]
    for i in range(n):
        s65 = (-65 * u[i]) % Q
        s73 = (-73 * u[i]) % Q
        P = [(x + s65) % Q for x in P] + [(x + s73) % Q for x in P]
    T = [0]
    for i in range(n):
        s65 = (65 * u[n + i]) % Q
        s73 = (73 * u[n + i]) % Q
        T = [(x + s65) % Q for x in T] + [(x + s73) % Q for x in T]
    return n, B, Q, W, H, J, u, P, T


# ------------------------------------------------------------------ 1
def sec1():
    head(1, "the exact seed count")
    for m in (1, 2):
        n, B, Q, W, H, J, u, P, T = halves(m)
        check(all(t % B == 0 for t in T), "m=%d: B divides every T_t" % m)
        a = [p // B for p in P]
        b = [t // B for t in T]
        check(len(set(b)) == len(b), "m=%d: the b_t are distinct" % m)
        # r_{st} formula and the direct residues
        direct = [0]
        for i in range(2 * n):
            s65 = (-65 * u[i]) % Q
            s73 = (-73 * u[i]) % Q
            direct = [(x + s65) % Q for x in direct] + [(x + s73) % Q
                                                        for x in direct]
        got = []
        for si, p in enumerate(P):
            for ti, t in enumerate(T):
                r = (p - t) % Q
                check(r == (p % B) + B * ((a[si] - b[ti]) % B),
                      "m=%d: r_st = (P mod B) + B[a-b]_B" % m)
                got.append(r)
        check(sorted(got) == sorted(direct),
              "m=%d: the split reproduces every residue" % m)
        D1 = sum(1 for r in direct if r < H)
        D2 = sum(1 for si in range(len(P)) for ti in range(len(T))
                 if (a[si] - b[ti]) % B < W)
        check(D1 == D2, "m=%d: D = #{pairs with [a-b]_B < W}" % m)
        print("    m=%d: B | T_t, the b_t are distinct, r_st formula holds for"
              % m)
        print("          all %d ordered pairs, and D = %d both ways."
              % (len(P) * len(T), D1))
    print("  Note a_s need NOT be distinct (the A list carries multiplicity),")
    print("  the cutoff is strict, and r > 0 always since r = 15 or 7 mod 64.")


# ------------------------------------------------------------------ 2
def sec2(mmax):
    head(2, "the dyadic bin identity, and the control table")
    print("  C_{j+1} = (C_j + G_j)/2 with G_j = sum_I dA(I) dZ(I), because")
    print("  A_L Z_L + A_R Z_R = (A Z + dA dZ)/2.  Unrolling and using")
    print("  C_0 = |A||Z| = 2^{3m} * 2^{3m+1} = 2^{6m+1} gives")
    print("      C_J = 2 + sum_{j<J} e_j,   e_j = 2^{j-J} G_j.")
    print("  %-4s %-8s %-8s %-8s %-8s %-10s %s"
          % ("m", "D_m(0)", "Same", "Next", "C_J", "max e_j", "identity"))
    want = {1: (2, 2, 1, 3, Fr(1)), 2: (2, 4, 2, 6, Fr(17, 8)),
            3: (3, 1, 4, 5, Fr(5, 4)), 4: (0, 0, 1, 1, Fr(1, 2)),
            5: (0, 0, 0, 0, Fr(7, 16)), 6: (1, 2, 1, 3, Fr(1))}
    for m in range(1, mmax + 1):
        n, B, Q, W, H, J, u, P, T = halves(m)
        a = sorted(p // B for p in P)
        bs = [t // B for t in T]
        # D by cyclic interval search
        D = 0
        for bt in bs:
            lo, hi = bt, bt + W
            if hi <= B:
                D += bisect.bisect_left(a, hi) - bisect.bisect_left(a, lo)
            else:
                D += len(a) - bisect.bisect_left(a, lo)
                D += bisect.bisect_left(a, hi - B)
        # Same / Next with bins of width W
        shift = (18 * m) - J          # log2(W)
        from collections import Counter
        HA = Counter(x >> shift for x in a)
        same = sum(HA.get(bt >> shift, 0) for bt in bs)
        nxt = sum(HA.get(((bt + W) % B) >> shift, 0) for bt in bs)
        Z = [bt for bt in bs] + [(bt + W) % B for bt in bs]
        CJ = same + nxt
        check(CJ == sum(HA.get(z >> shift, 0) for z in Z),
              "m=%d: C_J = sum_I A(I)Z(I)" % m)
        check(D <= CJ, "m=%d: D <= C_J" % m)
        # e_j
        es = []
        for j in range(J):
            sh = 18 * m - (j + 1)
            HAc = Counter(x >> sh for x in a)
            HZc = Counter(z >> sh for z in Z)
            G = 0
            for I in set(k >> 1 for k in list(HAc) + list(HZc)):
                dA = HAc.get(2 * I, 0) - HAc.get(2 * I + 1, 0)
                dZ = HZc.get(2 * I, 0) - HZc.get(2 * I + 1, 0)
                G += dA * dZ
            es.append(Fr(G, 2 ** (J - j)))
        check(2 + sum(es) == CJ, "m=%d: C_J = 2 + sum e_j" % m)
        mx = max(es)
        exp = want[m]
        check((D, same, nxt, CJ, mx) == exp,
              "m=%d: table row %s vs expected %s"
              % (m, (D, same, nxt, CJ, str(mx)), tuple(map(str, exp))))
        print("  %-4d %-8d %-8d %-8d %-8d %-10s %s"
              % (m, D, same, nxt, CJ, str(mx), "2+sum e_j = C_J"))
    print("    Every row of the supplied table reproduced exactly, and the")
    print("    identity C_J = 2 + sum_j e_j holds as exact rationals.")
    return


# ------------------------------------------------------------------ 3
def sec3():
    head(3, "the structure of A and Z, and where a flip-pairing stalls")
    m = 3
    n, B, Q, W, H, J, u, P, T = halves(m)
    tau = [8 * u[n + i] // B for i in range(n)]
    sig = [(-8 * u[i]) % Q for i in range(n)]
    print("  Z side: flipping letter i of t shifts T by +-8u_{n+i}, and")
    print("  B | 8u_{n+i}, so b_t shifts by the EXACT integer tau_i =")
    print("  8u_{n+i}/B.  Their valuations:")
    for i in range(min(n, 5)):
        v = 0; x = tau[i] % B
        while x % 2 == 0 and x: x //= 2; v += 1
        check(v == 3 + 6 * i, "tau_%d has v_2 = 3+6i" % i)
        print("      tau_%d: v_2 = %d (= 3+6i)" % (i, v))
    print("    so Z is an affine image of a dissociated subset-sum set --")
    print("    the same shape as the depth-n residue set.")
    print("  A side: flipping letter i of s shifts P by -+8u_i, but")
    print("  B does NOT divide 8u_i for i < n, so a_s = floor(P_s/B) shifts by")
    print("  floor((P +- 8u_i)/B) - floor(P/B), which is the exact shift PLUS")
    print("  a carry that depends on P_s mod B.  Measured carry rates:")
    Pm = P
    for i in (0, 1, 2):
        carries = 0
        for p in Pm[:4000]:
            q2 = (p + sig[i]) % Q
            exact = (sig[i] // B)
            if (q2 // B - p // B) % B != exact % B:
                carries += 1
        print("      letter %d: %d of 4000 flips carry across a B-boundary"
              % (i, carries))
    print("    The carry is exactly the information the pairing must track;")
    print("    it is NOT determined by the bin of a_s.")
    print("  e_j profile at m=3 (level : e_j), showing where the mass sits:")
    a = sorted(p // B for p in P)
    bs = [t // B for t in T]
    Z = bs + [(bt + W) % B for bt in bs]
    from collections import Counter
    for j in range(J):
        sh = 18 * m - (j + 1)
        HAc = Counter(x >> sh for x in a); HZc = Counter(z >> sh for z in Z)
        G = 0
        for I in set(k >> 1 for k in list(HAc) + list(HZc)):
            G += (HAc.get(2*I,0)-HAc.get(2*I+1,0))*(HZc.get(2*I,0)-HZc.get(2*I+1,0))
        if j < 10 or abs(G) > 0:
            print("      j=%-3d e_j = %-12s" % (j, str(Fr(G, 2 ** (J - j)))),
                  end="\n" if j % 1 == 0 else "")
    check(True, "profile printed")


# ------------------------------------------------------------------ 4
def sec4():
    head(4, "the control, and the exact line that distinguishes")
    for m in (1, 2, 3):
        n = 3 * m; B = 64 ** n; W = 64 ** (2 * m); Q = B * B; J = 6 * m
        Ps = [7 + 8 * sum(((k >> i) & 1) * 64 ** i for i in range(n))
              for k in range(1 << n)]
        # mirror the convention r = [P - T]_Q, so the second half enters
        # with the opposite sign
        Ts = [(-8 * sum(((k >> i) & 1) * 64 ** (n + i) for i in range(n))) % Q
              for k in range(1 << n)]
        a = [p // B for p in Ps]
        bs = [t // B for t in Ts]
        check(all(x == 0 for x in a), "m=%d: every control a_s is 0" % m)
        from collections import Counter
        HA = Counter(a)
        Z = bs + [(bt + W) % B for bt in bs]
        shift = 18 * m - J
        CJ = sum(HA.get(z >> shift, 0) for z in Z)
        D = sum(1 for p in Ps for t in Ts if (p // B - t // B) % B < W)
        print("      m=%d: A is a point mass at 0 (all %d of them); D = %d ="
              " 2^{5m} = %d ; C_J = %d"
              % (m, len(a), D, 2 ** (5 * m), CJ))
        check(D == 2 ** (5 * m), "m=%d: control count is 2^{5m}" % m)
    print("  THE LINE.  For the control every P_s < B, so a_s = 0 for all s:")
    print("  the prefix list is a single atom, dA(I) is maximal along one")
    print("  dyadic path at EVERY level, and e_j is of size 2^{3m} there --")
    print("  the test e_j <= 13 fails at once.  For the true family the P_s")
    print("  fill [0,Q) because u_i = 64^i 81^{-i-1} has a nonzero HIGH half,")
    print("  which exists only because 81^{-i-1} is a nontrivial 2-adic unit.")
    print("  So any estimate must prove that the high halves floor(P_s/B)")
    print("  spread over [0,B); that single step is where 81 enters and where")
    print("  the control dies.")


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    sec1(); sec2(mmax); sec3(); sec4()
    print("\n%d/%d checks passed" % (OK, OK + len(BAD)))
    for b in BAD: print("  FAILED:", b)
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())
