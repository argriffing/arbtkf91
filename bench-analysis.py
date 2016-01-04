from __future__ import print_function, division

import os
from subprocess import Popen, PIPE
import json

import data_source


def bench_pair(bench, model_params, a, b):
    d = model_params.copy()
    d.update(precision='float', samples=10, sequence_a=a, sequence_b=b)
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
    """
    print('verified:', out['verified'])
    if not out['verified']:
        print('input and output corresponding to verification failure:')
        print(s_in)
        print(stdout_data)
        raise Exception('verification failed')
    """
    #print('json stuff:')
    #print(out)
    elapsed = out['elapsed_ticks']
    print(elapsed)
    return sum(elapsed)


def main():
    bench = os.path.realpath('./bin/arbtkf91-bench')
    model_params = dict(
            pa_n=27, pa_d=100,
            pc_n=24, pc_d=100,
            pg_n=26, pg_d=100,
            pt_n=23, pt_d=100,
            lambda_n=1, lambda_d=1,
            mu_n=2, mu_d=1,
            tau_n=1, tau_d=10)
    for name, fin in data_source.gen_files():
        print(name)
        total_ticks = 0
        for a, b in data_source.gen_sequence_pairs(fin):
            total_ticks += bench_pair(bench, model_params, a, b)
        print(total_ticks / 10)
        print()


if __name__ == '__main__':
    main()
