"""Energies benchmark entry point: python -m test.benchmarks.energies <tier> [options].

T1 (regression) and T2 to T5 (accuracy) run from here (see README.md).
"""
import sys
from functools import partial
from test.benchmarks.energies import t1, t2

TIERS = {'t1': t1.main, 't2': partial(t2.main, tier='t2'), 't3': partial(t2.main, tier='t3'),
         't4': partial(t2.main, tier='t4'), 't5': partial(t2.main, tier='t5')}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in TIERS:
        print(__doc__)
        print('Tiers:', ', '.join(TIERS))
        sys.exit(2)
    sys.exit(TIERS[sys.argv[1]](sys.argv[2:]))
