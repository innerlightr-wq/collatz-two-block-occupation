import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prefix6m import starts_below
m = 7
h = starts_below(6 * m, 64 ** (5 * m))
print("m=%d (t=%d): matches = %d, allowance 84m = %d" % (m, 6 * m, len(h), 84 * m))
for x in h:
    print("  ", x)
