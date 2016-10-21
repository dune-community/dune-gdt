#!/usr/bin/env python3

import os
import pickle
import sys
from pprint import pprint
import subprocess
import time
from contextlib import contextmanager
import binpacking
from multiprocessing import Pool


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


def _compile(binary):
    with elapsed_timer() as timer:
        subprocess.check_call(['ninja', '-j1', binary])
        return timer()


def _run_tests(tpl):
    binary, teststrings = tpl
    testtimes = 0
    for test in teststrings.split('/'):
        with elapsed_timer() as timer:
            subprocess.check_call(['ctest', '-j1', '-R', test])
            testtimes += timer()
        return testtimes


def redo_timings(builddir, binaries, testnames, processes):
    os.chdir(builddir)
    testlimit = -1
    binaries = binaries[:testlimit]
    testnames = testnames[:testlimit]

    with Pool(processes=processes) as pool:
        compiles = pool.map(_compile, binaries)
    with Pool(processes=processes) as pool:
        testruns = pool.map(_run_tests, zip(binaries, testnames))

    totals = [a+b for a,b in zip(compiles, testruns)]

    print('totals')
    pprint(totals)
    pickle.dump(totals, open(os.path.join(builddir, pickle_file), 'wb'))
    return {b: t for b, t in zip(binaries, totals)}


# list comes with a leading empty entry
testnames = sys.argv[4].split('/')[1:]
builddir = sys.argv[1]
binaries = sys.argv[3].split(';')
processes = 4

try:
    totals = pickle.load(open(os.path.join(builddir, pickle_file), 'rb'))
    if totals.keys() != binaries:
        totals = redo_timings(builddir, binaries, testnames, processes)
except FileNotFoundError:
    totals = redo_timings(builddir, binaries, testnames, processes)

builder_count = sys.argv[2]

b = list(totals.keys())
bins = binpacking.to_constant_volume(b,MAXTIME)
pprint(bins)