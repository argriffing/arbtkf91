from __future__ import print_function, division

import argparse
import os
from subprocess import Popen, PIPE
import json

import data_source

def rat(a, b):
    return dict(num=a, denom=b)

def bench_pair(bench, d):
    s_in = json.dumps(d)
    args = [bench]
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    out = json.loads(stdout_data)
    return out['elapsed_ticks']


def main(args):
    bench = 'arbtkf91-bench'
    data_path = os.path.realpath(args.bench_data)
    samples = args.samples
    model_params = {
            "pa" : rat(27, 100),
            "pc" : rat(24, 100),
            "pg" : rat(26, 100),
            "pt" : rat(23, 100),
            "lambda" : rat(1, 1),
            "mu" : rat(2, 1),
            "tau" : rat(1, 10)}
    for name, fin in data_source.gen_files(data_path):
        print(name)
        total_ticks = 0
        for a, b in data_source.gen_sequence_pairs(fin):
            d = dict(
                precision=args.precision,
                rtol=0.0,
                samples=samples,
                parameters=model_params,
                sequence_a=a,
                sequence_b=b)
            elapsed = bench_pair(bench, d)
            print(len(a), len(b), elapsed)
            total_ticks += sum(elapsed)
        print(total_ticks / samples)
        print()


def pos_int(x):
    x = int(x)
    assert x > 0
    return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bench-data', required=True,
            help='benchmark data directory')
    parser.add_argument('--samples', default=10, type=pos_int,
            help='number of runs per alignment for benchmarking')
    parser.add_argument('--precision', required=True,
            choices=('float', 'double', 'mag', 'arb256'))
    args = parser.parse_args()
    main(args)
