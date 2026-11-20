#!/usr/bin/env python3
"""Independent confirmation that r_min(49) is minimal.

Instead of a successor search for the minimum, this counts EVERY 49-block
start below the claimed minimum plus one, by an exact cyclic-interval test,
and with the COMPLEMENTARY split of the letters (store 0..23, iterate
24..48).  The answer must be the single claimed value.

Resources: ~1.3 GiB resident, a few minutes.
"""
import bisect
import sys
import time

Q = 49
CLAIM = 165388305886088784630384508139202367487215600393310936214638366216510132943
C = {"A": 65, "B": 73}


def main():
    t0 = time.time()
    M = 64 ** Q
    mask = M - 1
    u = [(64 ** i * pow(81, -i - 1, M)) & mask for i in range(Q)]
    base = (-65 * sum(u)) & mask
    s = [(-8 * x) & mask for x in u]

    def sums(idx):
        out = [0]
        for i in idx:
            out = out + [(a + s[i]) & mask for a in out]
        return out

    store = sums(range(0, 24))
    store.sort()
    n = len(store)
    print("stored %d subset sums of letters 0..23 in %.1fs"
          % (n, time.time() - t0))
    sys.stdout.flush()
    A1 = sums(range(24, 36))
    A2 = sums(range(36, Q))
    limit = CLAIM + 1
    print("iterating %d x %d subset sums of letters 24..48, cutoff = claim + 1"
          % (len(A1), len(A2)))
    sys.stdout.flush()

    hits = []
    bl = bisect.bisect_left
    for a1 in A1:
        v1 = a1
        for a2 in A2:
            lo = (-base - v1 - a2) & mask
            hi = lo + limit
            if hi <= M:
                j = bl(store, lo)
                while j < n and store[j] < hi:
                    hits.append((base + v1 + a2 + store[j]) & mask)
                    j += 1
            else:
                j = bl(store, lo)
                while j < n:
                    hits.append((base + v1 + a2 + store[j]) & mask)
                    j += 1
                j = 0
                hi -= M
                while j < n and store[j] < hi:
                    hits.append((base + v1 + a2 + store[j]) & mask)
                    j += 1
    hits.sort()
    print("49-block starts below claim+1: %d" % len(hits))
    for h in hits:
        print("   %d" % h)
    ok = hits == [CLAIM]
    print("exactly the claimed minimum: %s" % ok)
    # replay for completeness
    x, word = CLAIM, ""
    for _ in range(Q):
        r = x % 64
        lab = "A" if r == 15 else ("B" if r == 7 else None)
        assert lab is not None
        word += lab
        x = (81 * x + C[lab]) // 64
    print("replay: %d blocks, word %s, endpoint mod 64 = %d"
          % (len(word), word, x % 64))
    print("elapsed %.1fs" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
