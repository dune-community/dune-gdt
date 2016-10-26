#!/usr/bin/env python3

import os
import pickle
import sys
from pprint import pprint
import subprocess
import time
from contextlib import contextmanager
import binpacking
from multiprocessing import Pool, cpu_count


MAXTIME = 4*60
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
        try:
            _ = subprocess.check_output(['ninja', '-j1', binary], universal_newlines=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as cpe:
            if 'Timeout' not in cpe.output:
                raise cpe
            print('Timeout in compile {}'.format(binary))
        return timer()


def _run_tests(tpl):
    binary, teststrings = tpl
    testtimes = 0
    for test in teststrings.split(';'):
        with elapsed_timer() as timer:
            try:
                _ = subprocess.check_output(['ctest', '-j1', '-N', '-R', test], universal_newlines=True,
                                            stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as cpe:
                if 'Timeout' not in cpe.output:
                    raise cpe
                # be pessimistic and double the timeout value as time for this run
                testtimes += timer()
                print('Timeout in {} from {}'.format(test, binary))
            testtimes += timer()
    return testtimes


def _redo(processes, keys, *args):
    try:
        with Pool(processes=processes) as pool:
            result = pool.map(*args)
        return {k: v for k,v in zip(keys, result)}
    except subprocess.CalledProcessError as cpe:
        print('*'*79)
        print(cpe.stdout)
        print(cpe.stderr)
        print('*' * 79)
        raise cpe

def do_timings(builddir, pickledir, binaries, testnames, processes):
    os.chdir(builddir)
    testlimit = -1

    binaries = binaries[:testlimit]
    compiles_fn = os.path.join(pickledir, 'compiles_' + pickle_file)
    try:
        compiles = pickle.load(open(compiles_fn, 'rb'))
        if set(compiles.keys()) != set(binaries):
            print('redoing compiles due to mismatched binaries')
            compiles = _redo(processes, binaries, _compile, binaries)
    except FileNotFoundError:
        print('redoing compiles due to missing pickle')
        compiles = _redo(processes, binaries, _compile, binaries)
    pickle.dump(compiles, open(compiles_fn, 'wb'))

    testnames = testnames[:testlimit]
    testruns_fn = os.path.join(pickledir, 'testruns_' + pickle_file)
    try:
        loaded_testnames, testruns = pickle.load(open(testruns_fn, 'rb'))
        if set(compiles.keys()) != set(binaries) or loaded_testnames != testnames:
            print('redoing tests due to mismatched binaries/testnames')
            testruns = _redo(processes, binaries, _run_tests, zip(binaries, testnames))
    except FileNotFoundError:
        print('redoing tests due to missing pickle')
        testruns = _redo(processes, binaries, _run_tests, zip(binaries, testnames))
    pickle.dump((testnames, testruns), open(testruns_fn, 'wb'))

    totals = {n: compiles[n]+testruns[n] for n in binaries}
    pickle.dump(totals, open(os.path.join(pickledir, pickle_file), 'wb'))
    print('totals')
    pprint(totals)
    return totals


# list comes with a leading empty entry
testnames = sys.argv[4].split('/')[1:]
builddir = sys.argv[1]
pickledir = sys.argv[2]
binaries = sys.argv[3].split(';')
processes = cpu_count()

totals = do_timings(builddir, pickledir, binaries, testnames, processes)
builder_count = sys.argv[2]

b = list(totals.keys())
bins = binpacking.to_constant_volume(totals, MAXTIME)
for idx, bin in enumerate(bins):
    pprint('Bin {} vol: {}'.format(idx, sum(bin.values())))
    pprint(bin)