#!/usr/bin/env python3
"""Exact least positive seed completing 49 aligned two-block steps.

r_min(q) = min over w in {A,B}^q of the canonical representative of
  r_w = -81^{-q} C_w  (mod 64^q),    C_w = sum_i 81^{q-1-i} 64^i c_{w_i}.
The 2^q residues are base + sum_{i in B(w)} s_i mod M, so the minimum is
found by a cyclic successor search: store the subset sums of one half,
sorted, and for each subset sum a of the other half take the least stored
b >= M - (base+a) (the wrapped candidate) and the least stored b overall
(the unwrapped candidate).  Every word is one (a, b) pair, and for each a
the b minimising (a+b) mod M is among those two, so the search is
exhaustive over all 2^q words.

Split used here: STORE letters 25..48, ITERATE letters 0..24.
Resources: ~1.3 GiB resident, a few minutes.
"""
import bisect
import sys
import time

Q = 49
EXPECTED = 165388305886088784630384508139202367487215600393310936214638366216510132943
C = {"A": 65, "B": 73}


def main():
    t0 = time.time()
    M = 64 ** Q
    mask = M - 1
    u = [(64 ** i * pow(81, -i - 1, M)) & mask for i in range(Q)]
    base = (-65 * sum(u)) & mask
    s = [(-8 * x) & mask for x in u]
    print("M = 64^%d = 2^%d" % (Q, 6 * Q))

    def sums(idx):
        out = [0]
        for i in idx:
            out = out + [(a + s[i]) & mask for a in out]
        return out

    store = sums(range(25, Q))
    store.sort()
    n = len(store)
    b0 = store[0]
    print("stored %d subset sums of letters 25..48 in %.1fs"
          % (n, time.time() - t0))
    sys.stdout.flush()

    A1 = sums(range(0, 12))
    A2 = sums(range(12, 25))
    print("iterating %d x %d = %d subset sums of letters 0..24"
          % (len(A1), len(A2), len(A1) * len(A2)))
    sys.stdout.flush()

    best = M
    bl = bisect.bisect_left
    done = 0
    for a1 in A1:
        v1 = (base + a1) & mask
        for a2 in A2:
            v = (v1 + a2) & mask
            # unwrapped candidate
            c = v + b0
            if c < M:
                if c < best:
                    best = c
            # wrapped candidate: least stored b >= M - v
            j = bl(store, M - v)
            if j < n:
                c = v + store[j] - M
                if c < best:
                    best = c
        done += len(A2)
        if done % (1 << 22) == 0:
            print("  %d/%d done, best so far %d (%.1fs)"
                  % (done, len(A1) * len(A2), best, time.time() - t0))
            sys.stdout.flush()

    print("\nr_min(%d) = %d" % (Q, best))
    print("digits: %d" % len(str(best)))
    print("expected match: %s" % (best == EXPECTED))
    # replay
    x, word = best, ""
    for _ in range(Q):
        r = x % 64
        lab = "A" if r == 15 else ("B" if r == 7 else None)
        if lab is None:
            print("REPLAY FAILED at block %d: residue %d" % (len(word), r))
            return 1
        word += lab
        x = (81 * x + C[lab]) // 64
    print("word   : %s" % word)
    print("endpoint %d" % x)
    print("endpoint mod 64 = %d  (exits: %s)" % (x % 64, x % 64 not in (7, 15)))
    print("elapsed %.1fs" % (time.time() - t0))
    return 0 if best == EXPECTED else 1


if __name__ == "__main__":
    sys.exit(main())
