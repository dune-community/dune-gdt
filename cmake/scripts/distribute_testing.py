#!/usr/bin/env python3

import os
import pickle
import sys
from pprint import pprint
import subprocess
import time
from contextlib import contextmanager
import binpacking

MAXTIME = 45*60
pickle_file = 'totals.pickle'


@contextmanager
def elapsed_timer():
    clock = time.time
    start = clock()
    elapser = lambda: clock() - start
    yield lambda: elapser()
    end = clock()
    elapser = lambda: end-start


def redo_timings(argv):
    builddir = argv[1]
    builder_count = argv[2]
    binaries = argv[3].split(';')
    # list comes with a leading empty entry
    testnames = argv[4].split('/')[1:]

    os.chdir(builddir)
    compiles = {}
    testtimes = {}

    testlimit = -1
    binaries = binaries[:testlimit]
    testnames = testnames[:testlimit]
    totals = {}

    for binary in binaries:
        with elapsed_timer() as timer:
            subprocess.check_call(['ninja', '-j1', binary])
            compiles[binary] = timer()
            totals[binary] = compiles[binary]
    for binary, teststrings in zip(binaries, testnames):
        testtimes[binary] = 0
        pprint(teststrings.split('/'))
        for test in teststrings.split('/'):
            with elapsed_timer() as timer:
                subprocess.check_call(['ctest', '-j1', '-R', test])
                testtimes[binary] += timer()
                totals[binary] = testtimes[binary]
    print('compiles')
    pprint(compiles)
    print('testtimes')
    pprint(testtimes)
    print('totals')
    pprint(totals)
    pickle.dump(totals, open(os.path.join(builddir, pickle_file), 'wb'))
    return totals


try:
    builddir = sys.argv[1]
    totals = pickle.load(open(os.path.join(builddir, pickle_file), 'rb'))
except FileNotFoundError:
    totals = redo_timings(sys.argv)

b = list(totals.values())
bins = binpacking.to_constant_volume(b,MAXTIME)

pprint(bins)