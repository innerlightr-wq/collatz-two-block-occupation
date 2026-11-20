#!/usr/bin/env python3
"""Checks for the claims new in revision 2 of the occupation note."""
import math, sys
from collections import Counter
from fractions import Fraction as Fr
OK, BAD = 0, []
def ck(c, m):
    global OK
    if c: OK += 1
    else: BAD.append(m); print("FAIL:", m)

# --- Remark 2.5: the floor retains the additive constant of Theorem 2.2
RHO = Fr(81, 64)
for t, b in ((12, 6), (40, 1), (48, 1)):
    H = 2 * 32 ** t
    exact = (Fr(H) + Fr(73, 17)) * RHO ** (-b) - Fr(73, 17)
    ck(exact > 0, "t=%d: the Theorem 2.2 floor is positive" % t)
lam32 = 5 * math.log(2) / math.log(81 / 64)
ck(abs(lam32 - 14.712372983) < 1e-8, "log_rho 32 = 14.712372983")
ck(3 ** 56 < 2 ** 89, "3^56 < 2^89 certifies rho^14 < 32")
ck(abs(float(Fr(2) ** 89 / Fr(3) ** 56) - 1.1827129437) < 1e-9,
   "2^89/3^56 = 1.1827129437")
# the Remark-2.5 floor with its constant, at a sample depth
t = 100
floor_with_const = 2 * (Fr(2) ** 89 / Fr(3) ** 56) ** t - Fr(73, 17)
ck(floor_with_const > 0, "the exponential floor minus 73/17 is positive at t=100")

# --- moduli
ck(64 ** 48 == 2 ** 288 and 64 ** 49 == 2 ** 294, "64^48 = 2^288, 64^49 = 2^294")

# --- 5 log_rho 64
c88 = 5 * math.log(64) / math.log(81 / 64)
ck(abs(c88 - 88.274237898) < 1e-7, "5 log_rho 64 = 88.274237898")

# --- Theorem (depth averaging): the four ingredients
for j in (12, 16, 24):
    o, x = 1, 81 % (1 << j)
    while x != 1:
        x = x * 81 % (1 << j); o += 1
    ck(o == 2 ** (j - 4), "order of 81 mod 2^%d is 2^{j-4}" % j)
for (l, s) in ((1, 1), (2, 2)):
    M = 64 ** (l + s); L = 64 ** l
    ck(2 ** (2 * (l + s) - 1) * 6 * L == 3 * 2 ** (2 * l - 4 * s) * M,
       "l=%d s=%d: union bound equals 3*2^{2l-4s}" % (l, s))
    ck(2 ** (6 * l - 3 + 4) <= 2 * L, "l=%d: 2^{v+4} <= 2L at the top v" % l)
ck(3 * 2.0 ** (-2 * 1) == 0.75 and abs(3 * 2.0 ** (-4) - 0.1875) < 1e-12,
   "the stated bad-proportion bounds 0.75 and 0.1875")

# --- Proposition (half-word reduction) and the worked example at m=2, j=9
m = 2; n = 3 * m; B = 64 ** n; Q = B * B; W = 64 ** (2 * m); J = 6 * m
u = [(64 ** i * pow(81, -i - 1, Q)) % Q for i in range(2 * n)]
P = [0]
for i in range(n):
    P = [(x - 65 * u[i]) % Q for x in P] + [(x - 73 * u[i]) % Q for x in P]
T = [0]
for i in range(n):
    T = [(x + 65 * u[n + i]) % Q for x in T] + [(x + 73 * u[n + i]) % Q for x in T]
ck(all(t % B == 0 for t in T), "B divides every T_t")
a = sorted(p // B for p in P); bs = [t // B for t in T]
Z = bs + [(b + W) % B for b in bs]
ck(len(a) == 2 ** n and len(Z) == 2 ** (n + 1), "|A| = 2^{3m}, |Z| = 2^{3m+1}")
D = sum(1 for p in P for t in T if (p // B - t // B) % B < W)
ck(D == 2, "m=2: D_m(0) = 2 by the half-word count")
tot = Fr(0); e9 = None; sdz9 = None
for j in range(J):
    sh = 6 * n - (j + 1)
    HA = Counter(x >> sh for x in a); HZ = Counter(z >> sh for z in Z)
    par = set(k >> 1 for k in list(HA) + list(HZ))
    G = sum((HA.get(2 * I, 0) - HA.get(2 * I + 1, 0)) *
            (HZ.get(2 * I, 0) - HZ.get(2 * I + 1, 0)) for I in par)
    tot += Fr(G, 2 ** (J - j))
    if j == 9:
        e9 = Fr(G, 2 ** (J - j))
        sdz9 = sum(HZ.get(2 * I, 0) - HZ.get(2 * I + 1, 0) for I in par)
        nzZ = sum(1 for I in par if HZ.get(2*I,0) != HZ.get(2*I+1,0))
        npar = len(par)
ck(2 + tot == 6, "m=2: C_J = 2 + sum e_j = 6")
ck(sdz9 == 0, "m=2, j=9: sum_I dZ(I) = 0")
ck(e9 == Fr(17, 8), "m=2, j=9: e_9 = 17/8")
ck(npar == 110 and nzZ == 54, "m=2, j=9: 110 parents, 54 with dZ != 0")
ck(B == 2 ** 36 and W == 2 ** 24, "m=2: B = 2^36, W = 2^24")
# tau_i exactness
tau = [8 * u[n + i] // B for i in range(n)]
for i in range(n):
    x = tau[i]; v = 0
    while x % 2 == 0: x //= 2; v += 1
    ck(v == 3 + 6 * i, "v_2(tau_%d) = 3+6i" % i)
    ck(8 * u[n + i] % B == 0, "B divides 8u_{n+i} for i=%d" % i)
# --- the corrected pairing statement (m=1 counterexample)
m1 = 1; n1 = 3 * m1; B1 = 64 ** n1; Q1 = B1 * B1; W1 = 64 ** (2 * m1)
u1 = [(64 ** i * pow(81, -i - 1, Q1)) % Q1 for i in range(2 * n1)]
T1 = [0]
for i in range(n1):
    T1 = [(x + 65 * u1[n1 + i]) % Q1 for x in T1] + \
         [(x + 73 * u1[n1 + i]) % Q1 for x in T1]
b1 = [t // B1 for t in T1]
Z1 = set(b1 + [(b + W1) % B1 for b in b1])
tau0 = (8 * u1[n1] // B1) % B1
ck(B1 == 262144 and tau0 == 116232, "m=1: B = 262144, tau_0 = 116232")
ck(25921 in Z1 and (25921 + tau0) % B1 == 142153 and 142153 in Z1,
   "m=1: 25921 and 25921+tau_0 = 142153 both lie in Z")
ck((25921 + 2 * tau0) % B1 == 258385 and 258385 not in Z1,
   "m=1: 25921+2tau_0 = 258385 does NOT lie in Z (no translation invariance)")

# --- the corrected level counts at m=2, j=9
sh = 6 * n - 10
HA10 = Counter(x >> sh for x in a); HZ10 = Counter(z >> sh for z in Z)
children = set(list(HA10) + list(HZ10)); parents = set(k >> 1 for k in children)
ck(len(children) == 125, "m=2, j=9: 125 level-10 child bins hold an A or Z value")
ck(len(HZ10) == 78, "m=2, j=9: 78 level-10 child bins hold a Z value")
ck(len(parents) == 110, "m=2, j=9: 110 occupied level-9 parent bins")
ck(sum(1 for I in parents if HZ10.get(2*I,0) != HZ10.get(2*I+1,0)) == 54,
   "m=2, j=9: 54 parents with dZ != 0")

# --- aligned vs unaligned six-digit windows
import itertools as _it
_A, _B = "111100", "111010"
_counts = set()
for _w in _it.product((_A, _B), repeat=3):
    _s = "".join(_w)
    for _o in range(13):
        if len(_s[_o:_o+6]) == 6:
            _counts.add((_o % 6, _s[_o:_o+6].count("1")))
ck({c for o, c in _counts if o == 0} == {4},
   "every ALIGNED six-step block has exactly four ones")
ck({c for o, c in _counts} == {3, 4, 5},
   "unaligned windows take 3, 4 or 5 ones")
ck((_A + _B)[4:10] == "001110" and (_A + _B)[4:10].count("1") == 3,
   "001110 occurs in AB at offset 4 with three ones")

print("%d/%d checks passed" % (OK, OK + len(BAD)))
for b in BAD: print("  FAILED:", b)
sys.exit(1 if BAD else 0)
