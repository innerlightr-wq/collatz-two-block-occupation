#!/usr/bin/env python3
"""Every numerical claim carried into the manuscript, checked exactly.

Run: python3 note_checks.py [mmax]     (mmax defaults to 6; m=7,8 are cited
from m7.log and the depth-48 logs).  Needs ~1 GiB for mmax=6, no network.
"""
import bisect
import math
import sys
from fractions import Fraction as Fr

OK, BAD = 0, []
RHO = Fr(81, 64)
C = {"A": 65, "B": 73}
WORD_A, WORD_B = (1, 1, 1, 1, 0, 0), (1, 1, 1, 0, 1, 0)
MIN49 = 165388305886088784630384508139202367487215600393310936214638366216510132943
SEED48 = 100766337511639514688053192519038392361740472623921585501397359010728775


def check(cond, msg):
    global OK
    if cond:
        OK += 1
    else:
        BAD.append(msg)
        print("FAIL:", msg)


def head(s):
    print("\n--- " + s)


def parity_word(x, n=6):
    """n shortcut steps, parity recorded BEFORE each step."""
    bits = []
    for _ in range(n):
        p = x & 1
        bits.append(p)
        x = (3 * x + 1) // 2 if p else x // 2
    return tuple(bits), x


def label(x):
    r = x % 64
    return "A" if r == 15 else ("B" if r == 7 else None)


def step(x):
    return (81 * x + C[label(x)]) // 64


def switches(t):
    M = 64 ** t
    mask = M - 1
    u = [(64 ** i * pow(81, -i - 1, M)) & mask for i in range(t)]
    return (-65 * sum(u)) & mask, [(-8 * x) & mask for x in u], M


def residues(t):
    base, s, M = switches(t)
    mask = M - 1
    out = [base]
    for v in s:
        out = out + [(a + v) & mask for a in out]
    return out


def starts_below(t, limit):
    base, s, M = switches(t)
    mask = M - 1
    half = t // 2

    def subs(idx):
        out = [0]
        for i in idx:
            out = out + [(a + s[i]) & mask for a in out]
        return out
    SR = sorted(subs(range(half, t)))
    n = len(SR)
    hits = []
    for a in subs(range(half)):
        lo = (-base - a) & mask
        hi = (lo + limit) & mask
        if lo < hi:
            j = bisect.bisect_left(SR, lo)
            while j < n and SR[j] < hi:
                hits.append((base + a + SR[j]) & mask)
                j += 1
        else:
            j = bisect.bisect_left(SR, lo)
            while j < n:
                hits.append((base + a + SR[j]) & mask)
                j += 1
            j = 0
            while j < n and SR[j] < hi:
                hits.append((base + a + SR[j]) & mask)
                j += 1
    return sorted(hits)


# ---------------------------------------------------------------- section 1
def sec1():
    head("1. the two blocks, their classes, and the block maps")
    a = [r for r in range(64) if parity_word(r)[0] == WORD_A]
    b = [r for r in range(64) if parity_word(r)[0] == WORD_B]
    check(a == [15], "111100 occurs exactly on 15 mod 64 (got %s)" % a)
    check(b == [7], "111010 occurs exactly on 7 mod 64 (got %s)" % b)
    print("    A = 111100 <-> 15 mod 64 ;  B = 111010 <-> 7 mod 64")
    for k in range(0, 2000):
        x = 15 + 64 * k
        w, y = parity_word(x)
        check(w == WORD_A and y == (81 * x + 65) // 64 and (81 * x + 65) % 64 == 0,
              "A: x=%d gives (81x+65)/64" % x)
        x = 7 + 64 * k
        w, y = parity_word(x)
        check(w == WORD_B and y == (81 * x + 73) // 64 and (81 * x + 73) % 64 == 0,
              "B: x=%d gives (81x+73)/64" % x)
    print("    F_A(x) = (81x+65)/64 and F_B(x) = (81x+73)/64 verified on the")
    print("    first 2000 members of each class, by six-step simulation.")
    check(parity_word(15)[1] == 20 and parity_word(7)[1] == 10,
          "15 -> 20 and 7 -> 10")
    print("    15 -> 20 and 7 -> 10: the terminal parity is UNRESTRICTED, so")
    print("    completing a word says nothing about continuing.")
    # odd density
    check(sum(WORD_A) == 4 and sum(WORD_B) == 4, "each block has 4 odd steps")
    check(Fr(2, 3) > Fr(1, 1) * Fr(math.log(2) / math.log(3)).limit_denominator(10**9),
          "2/3 > log_3 2")
    print("    odd density = 4/6 = 2/3 = 0.666667 > log_3 2 = %.6f"
          % (math.log(2) / math.log(3)))
    # one class mod 64^q, by brute force
    for q in (1, 2, 3):
        res = residues(q)
        for n in range(1 << q):
            w = "".join("B" if (n >> i) & 1 else "A" for i in range(q))
            real = []
            for x in range(64 ** q):
                y, ok = x, True
                for j in range(q):
                    if label(y) != w[j]:
                        ok = False
                        break
                    y = step(y)
                if ok:
                    real.append(x)
            check(real == [res[n] % 64 ** q],
                  "q=%d word %s: realizers mod 64^q are exactly {r_w}"
                  % (q, w))
        print("    q=%d: all %d words -- realizer set is exactly one class"
              " mod 64^%d." % (q, 1 << q, q))


# ---------------------------------------------------------------- section 2
def sec2():
    head("2. the occupation inequality")
    for x in (7, 15, 71, 79, 2 ** 50 + 7, SEED48, MIN49):
        if label(x) is None:
            continue
        y = step(x)
        Y, Yp = 17 * x + 73, 17 * y + 73
        check(64 * Yp == 81 * Y - (136 if label(x) == "A" else 0),
              "64(17x'+73) = 81(17x+73) - %d at x=%d"
              % (136 if label(x) == "A" else 0, x))
        check(Fr(Yp) <= RHO * Y, "17x'+73 <= rho(17x+73) at x=%d" % x)
    print("    64(17x'+73) = 81(17x+73) for B and 81(17x+73)-136 for A:")
    print("    exact, so 17x_j+73 <= rho^j (17x_0+73) with equality iff all B.")
    # monotone r_min on the known values, and the floor formula
    for (t, b, H) in ((48, 1, 2 ** 241), (20, 0, 2 ** 101), (12, 6, 2 ** 61)):
        floor = (Fr(H) + Fr(73, 17)) * RHO ** (-b) - Fr(73, 17)
        print("    t=%-3d b=%-2d : r_min(t+b) >= %.6e" % (t, b, float(floor)))
        check(floor > 0, "t=%d: the floor is positive" % t)
    check(float((Fr(2 ** 241) + Fr(73, 17)) / RHO - Fr(73, 17)) > 2.79e72,
          "the depth-48 instance gives r_min(49) >= 2.79e72")
    # lambda and the 84 benchmark, exactly
    lam = 5 * math.log(2) / math.log(81 / 64)
    check(3 ** 56 < 2 ** 89, "rho^14 < 32 exactly (3^56 < 2^89)")
    check(3 ** 60 > 2 ** 95, "rho^15 > 32 exactly (3^60 > 2^95)")
    check(3 ** 336 < 2 ** 534, "rho^84 < 2^30 exactly (3^336 < 2^534)")
    print("    log_rho 32 = %.9f, 5 log_rho 64 = %.9f" % (lam, 6 * lam))
    print("    certified: 3^56 < 2^89 (rho^14 < 32), 3^336 < 2^534 (rho^84 <")
    print("    64^5), so 84 < 5 log_rho 64 = %.5f and 84m is admissible."
          % (6 * lam))
    check(abs(6 * lam - 88.274237898) < 1e-7, "5 log_rho 64 = 88.274237898")
    check(abs(lam - 14.712372983) < 1e-7, "log_rho 32 = 14.712372983")


# ---------------------------------------------------------------- section 3
def sec3():
    head("3. the depth-49 minimum and the depth-48 control")
    for (name, seed, q, endres) in (("r_min(49)", MIN49, 49, 31),
                                    ("the 2^241 seed", SEED48, 48, 35)):
        x, w = seed, ""
        while True:
            lab = label(x)
            if lab is None:
                break
            w += lab
            x = step(x)
        check(len(w) == q, "%s realizes exactly %d blocks" % (name, q))
        check(x % 64 == endres, "%s: endpoint is %d mod 64" % (name, endres))
        # replay through the raw shortcut steps
        y, labs = seed, []
        for _ in range(q):
            bits, y = parity_word(y)
            labs.append("A" if bits == WORD_A else ("B" if bits == WORD_B else "?"))
        check("".join(labs) == w, "%s: raw six-bit replay agrees" % name)
        check(y == x, "%s: raw replay reaches the same endpoint" % name)
        print("    %-15s %3d digits, %d blocks, endpoint %d mod 64, word"
              % (name, len(str(seed)), q, endres))
        print("        %s" % w)
    check(MIN49 > 2 ** 241, "r_min(49) exceeds 2^241, as the depth-48 "
                            "certificate requires")
    print("    consistency: r_min(49) = %.4e > 2^241 = %.4e."
          % (float(MIN49), 2.0 ** 241))


# ---------------------------------------------------------------- section 4
def sec4(mmax):
    head("4. the distinguished bin D_m(0), the mean-one lemma, and B_t")
    rows = []
    for m in range(1, mmax + 1):
        hits = starts_below(6 * m, 64 ** (5 * m))
        rows.append((m, len(hits)))
        check(len(hits) <= 84 * m, "m=%d: D_m(0) = %d <= 84m" % (m, len(hits)))
    print("    D_m(0) for m = 1..%d : %s"
          % (mmax, ", ".join(str(b) for _, b in rows)))
    check([b for _, b in rows][:6] == [2, 2, 3, 0, 0, 1],
          "D_m(0) = 2,2,3,0,0,1 for m = 1..6")
    print("    (m = 7 gives 2 and m = 8 gives 1; see m7.log and the depth-48")
    print("     certificate, where the unique start below 2^241 is < 2^240.)")
    # mean-one lemma over translates, exactly
    for m in (1, 2, 3):
        Hm = 64 ** m
        Sm = sorted(set(r % Hm for r in residues(m)))
        check(len(Sm) == 2 ** m, "m=%d: |S_m| = 2^m" % m)
        hist = {}
        for r in residues(5 * m):
            x, d = r, 0
            while d < 5 * m and label(x) is not None:
                x = step(x)
                d += 1
            check(d == 5 * m, "m=%d: every prefix representative realizes 5m"
                  % m)
            hist[x % Hm] = hist.get(x % Hm, 0) + 1
        tot = sum(hist.values())
        check(tot == 2 ** (5 * m), "m=%d: %d endpoints" % (m, tot))
        D = [sum(hist.get((c + t) % Hm, 0) for t in Sm) for c in range(Hm)]
        check(sum(D) == 2 ** (5 * m) * 2 ** m,
              "m=%d: sum over translates = 2^{5m} 2^m" % m)
        check(Fr(sum(D), Hm) == 1, "m=%d: the mean over translates is exactly 1"
              % m)
        d0 = [b for mm, b in rows if mm == m][0]
        check(D[0] == d0, "m=%d: D_m(0) from the translate table equals %d"
              % (m, d0))
        print("    m=%d: mean over all %d translates = 1 exactly; D_m(0) = %d,"
              " max_c D_m(c) = %d" % (m, Hm, D[0], max(D)))
    print("    so the average over bins is pinned, and the bin at the origin")
    print("    is one of 64^m translates: the mean cannot bound it.")
    # B_t census spot checks
    for t in (3, 12, 20):
        b = len(starts_below(t, 1 << (5 * t + 1)))
        check(b == {3: 2, 12: 6, 20: 0}[t], "B_%d = %d" % (t, b))
    print("    B_3 = 2, B_12 = 6, B_20 = 0 (full census t = 3..40 in"
          " windowcount.log: sum 65, max 6 at t = 12).")
    # comparison family
    for t in range(3, 15):
        K = -((-(5 * t - 2)) // 6)
        cnt = sum(1 for n in range(1 << t)
                  if 7 + 8 * sum(((n >> i) & 1) * 64 ** i for i in range(t))
                  < (1 << (5 * t + 1)))
        check(cnt == 2 ** K, "t=%d: comparison family count = 2^%d" % (t, K))
    print("    comparison family 7 + 8 sum eps_i 64^i: exactly")
    print("    2^ceil((5t-2)/6) members below 2^{5t+1}, for every t = 3..14.")


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    sec1()
    sec2()
    sec3()
    sec4(mmax)
    print("\n%d/%d checks passed" % (OK, OK + len(BAD)))
    for b in BAD:
        print("  FAILED:", b)
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())
