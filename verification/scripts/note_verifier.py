"""Exact verification of the 48-block count for the positive integer A/B model.

Run with Python 3.9 or later. No third-party libraries. The depth-48 run
sorts 2**24 Python integers and can use more than 1 GiB of memory.

For a word of t blocks, unrolling 64*x_next = 81*x + c gives the
starting residue modulo M = 64**t as sum(-c_i * 64**i / 81**(i+1)),
with the division performed modulo M. Each block congruence gives the
six prescribed shortcut input bits, without a terminal parity condition.

Start with an all-A word. Switching label i to B adds -8 * 64**i /
81**(i+1). With m=t//2, all late switches are multiples of K=64**m.
For an early sum a in [0,M), keep its low digit a mod K fixed. The late
switches add one of a sorted list of values b to the high digit a//K,
modulo M//K. The interval 0 <= x < 2**(5*t+1) is a union of complete
low-digit intervals. One modular interval lookup finds all late sums
which hit it. No rounding, sampling, or assumption about randomness.
"""

from bisect import bisect_left
from math import prod

A = (1, 1, 1, 1, 0, 0)
B = (1, 1, 1, 0, 1, 0)


def subset_sums(weights, mask):
    values = [0]
    for weight in weights:
        values += [(a + weight) & mask for a in values]
    return values


def count_starts(t):
    assert t >= 4 and t % 2 == 0
    midpoint = t // 2
    modulus = 64**t
    low_base = 64**midpoint
    high_base = modulus // low_base
    high_mask = high_base - 1
    full_mask = modulus - 1
    high_bits = high_mask.bit_length()
    low_bits = (low_base - 1).bit_length()
    limit = 1 << (5*t + 1)
    width, remainder = divmod(limit, low_base)
    assert remainder == 0 and width < high_base

    units = [(64**i * pow(81, -i-1, modulus)) & full_mask for i in range(t)]
    base = -65 * sum(units) & full_mask
    switches = [-8*u & full_mask for u in units]
    assert all(v % low_base == 0 for v in switches[midpoint:])
    right = subset_sums([v // low_base for v in switches[midpoint:]], high_mask)
    right.sort()
    size = len(right)

    # Bucket boundaries only accelerate bisect; comparisons are of full integers.
    bucket_bits = min(18, midpoint-1)
    shift = high_bits - bucket_bits
    buckets = [bisect_left(right, j << shift) for j in range(1 << bucket_bits)]
    buckets.append(size)

    part = midpoint // 2
    early1 = subset_sums(switches[:part], full_mask)
    early2 = subset_sums(switches[part:midpoint], full_mask)
    assert prod(map(len, [early1, early2, right])) == 2**t
    hits = []

    def record(a, b):
        x = (a + low_base*b) & full_mask
        assert 0 <= x < limit
        hits.append(x)

    for u in early1:
        base_u = (base + u) & full_mask
        for v in early2:
            a = (base_u + v) & full_mask
            start = -(a >> low_bits) & high_mask
            end = start + width
            bucket = start >> shift
            pos = bisect_left(right, start, buckets[bucket], buckets[bucket+1])
            while pos < size and right[pos] < end:
                record(a, right[pos])
                pos += 1
            if end > high_base:
                pos = 0
                stop = end - high_base
                while pos < size and right[pos] < stop:
                    record(a, right[pos])
                    pos += 1
    hits.sort()
    assert len(set(hits)) == len(hits)
    return hits


def shortcut_block(x):
    bits = []
    for _ in range(6):
        parity = x % 2
        bits.append(parity)
        x = (3*x + 1)//2 if parity else x//2
    return tuple(bits), x


def verify_replay(x, depth):
    for _ in range(depth):
        bits, x = shortcut_block(x)
        assert bits in (A, B)
    return x


def main():
    controls = {4: 560455, 10: 68794806996295, 12: 215776538488999239,
                20: 3281799071571577949650239917383}
    for t, expected in controls.items():
        # The depth-20 minimum exceeds 2**101. Check the empty interval there.
        hits = count_starts(t)
        if expected < 1 << (5*t+1):
            assert hits and hits[0] == expected
        else:
            assert hits == []
        for x in hits:
            verify_replay(x, t)
        print(f"Depth {t} passed. Starts in the test interval: {len(hits)}", flush=True)

    seed = 100766337511639514688053192519038392361740472623921585501397359010728775
    hits = count_starts(48)
    assert hits == [seed], hits
    endpoint = verify_replay(seed, 48)
    bits, _ = shortcut_block(endpoint)
    assert endpoint % 64 == 35 and bits not in (A, B)
    print(f"Only 48-block start in [1, 2**241): {seed}")
    print(f"Boundary after 48 blocks is {endpoint} (35 mod 64).")
    print("No positive seed below 2**241 completes 49 blocks.")


if __name__ == "__main__":
    main()
