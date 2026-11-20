#!/usr/bin/env python3
"""Independent count of 48-block starts below 2^241.

Different structure from the supplied verifier: branch on the B-positions
among the first 6 letters (64 branches), then meet-in-the-middle on two
blocks of 21 letters, with full 288-bit residues and an exact cyclic
interval test (no low/high digit decomposition).
"""
import bisect
import sys

T = 48
M = 64 ** T
MASK = M - 1
LIMIT = 1 << (5 * T + 1)          # 2^241


def units(t, mod):
    return [(64 ** i * pow(81, -i - 1, mod)) & (mod - 1) for i in range(t)]


def subset_sums(ws):
    out = [0]
    for w in ws:
        out = out + [(a + w) & MASK for a in out]
    return out


def count(t=T, nbranch=6, verbose=True):
    u = units(t, M)
    base = (-65 * sum(u)) & MASK
    s = [(-8 * x) & MASK for x in u]
    half = (t - nbranch) // 2
    midL = list(range(nbranch, nbranch + half))
    midR = list(range(nbranch + half, t))
    assert len(midL) + len(midR) + nbranch == t
    SR = subset_sums([s[i] for i in midR])
    SR.sort()
    n = len(SR)
    SL = subset_sums([s[i] for i in midL])
    hits = []
    for b in range(1 << nbranch):
        bb = base
        for i in range(nbranch):
            if b >> i & 1:
                bb = (bb + s[i]) & MASK
        for a in SL:
            lo = (-bb - a) & MASK
            hi = (lo + LIMIT) & MASK
            if lo < hi:
                p = bisect.bisect_left(SR, lo)
                while p < n and SR[p] < hi:
                    hits.append((bb + a + SR[p]) & MASK)
                    p += 1
            else:
                p = bisect.bisect_left(SR, lo)
                while p < n:
                    hits.append((bb + a + SR[p]) & MASK)
                    p += 1
                p = 0
                while p < n and SR[p] < hi:
                    hits.append((bb + a + SR[p]) & MASK)
                    p += 1
        if verbose and b % 8 == 7:
            print("  branch %d/%d done, hits so far %d" % (b + 1, 1 << nbranch,
                                                           len(hits)))
            sys.stdout.flush()
    return sorted(hits)


if __name__ == "__main__":
    hits = count()
    print("48-block starts below 2^241:", len(hits))
    for h in hits:
        print("  ", h)
    assert all(0 <= h < LIMIT for h in hits)
    print("EXPECTED one seed: 100766337511639514688053192519038392361740472623921585501397359010728775")
