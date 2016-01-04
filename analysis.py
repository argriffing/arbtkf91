from __future__ import print_function, division

import os
from subprocess import Popen, PIPE, STDOUT
import json

import numpy as np
from numpy.testing import assert_array_equal, assert_equal
from Bio import SeqIO

def select_indices(available, requested):
    incr = available // requested
    return [i*incr for i in range(requested)]


def analysis(bench, fin):
    fasta_objects = list(SeqIO.parse(fin, 'fasta'))
    sequences = [str(x.seq) for x in fasta_objects]
    available = len(sequences)
    requested = 10
    indices = select_indices(available, requested)
    selection = [sequences[i] for i in indices]
    assert_equal(len(selection), requested)

    model_params = dict(
            pa_n=27, pa_d=100,
            pc_n=24, pc_d=100,
            pg_n=26, pg_d=100,
            pt_n=23, pt_d=100,
            lambda_n=1, lambda_d=1,
            mu_n=2, mu_d=1,
            tau_n=1, tau_d=10)
    
    base = model_params.copy()
    base['precision'] = 'exact'
    base['samples'] = 10

    total_ticks = 0
    k = 0
    for i in range(requested):
        for j in range(i):
            a = selection[i]
            b = selection[j]

            d = base.copy()
            d['sequence_a'] = a
            d['sequence_b'] = b
            s_in = json.dumps(d)

            args = [bench]
            p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
            data = p.communicate(input=s_in)
            stdout_data, stderr_data = data
            try:
                out = json.loads(stdout_data)
            except ValueError as e:
                print('raw stdout data:')
                print(stdout_data)
                print('raw stderr data:')
                print(stderr_data)
                raise
            print('verified:', out['verified'])
            if not out['verified']:
                print('input and output corresponding to verification failure:')
                print(s_in)
                print(stdout_data)
                raise Exception('verification failed')
            #print('json stuff:')
            #print(out)
            k += 1
            elapsed = out['elapsed_ticks']
            print(elapsed)
            total_ticks += sum(elapsed)
    assert_equal(k, 45)
    print(total_ticks / 10)
    print()


def main():

    # mini test
    indices = select_indices(60, 10)
    assert_array_equal(indices[:3], [0, 6, 12])
    assert_array_equal(indices[-1:], [54])

    bench = os.path.realpath('./bin/arbtkf91-bench')
    check = os.path.realpath('./bin/arbtkf91-check')
    mypath = os.path.realpath('../../stamatakis/benchMark_data')
    for filename in os.listdir(mypath):
        if 'unaligned' in filename:
            print(filename)
            fullpath = os.path.join(mypath, filename)
            with open(fullpath) as fin:
                analysis(bench, fin)


if __name__ == '__main__':
    main()
