"""
Doing benchmark 0
Seq1 length: 10  Seq2 length: 10
Elapsed Time in millisecs: 0
Doing benchmark 1
Seq1 length: 100  Seq2 length: 100
Elapsed Time in millisecs: 0
Doing benchmark 2
Seq1 length: 1000  Seq2 length: 1000
Elapsed Time in millisecs: 5.5
Doing benchmark 3
Seq1 length: 10000  Seq2 length: 10000
Elapsed Time in millisecs: 771.6
Benchmarked succesfully!

"""
from __future__ import print_function, division

from subprocess import Popen, PIPE
import argparse
import os
import json

def rat(a, b):
    return dict(num=a, denom=b)

def bench_pair(bench, d):
    s_in = json.dumps(d)
    args = [bench]
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    out = json.loads(stdout_data)
    return 1000 * sum(out['elapsed_ticks']) / out['ticks_per_second']


def main(args):
    bench = 'arbtkf91-bench'
    with open(os.path.realpath(args.sequences)) as fin:
        lines = [x.strip() for x in fin.readlines()]
    npairs, remainder = divmod(len(lines), 2)
    if remainder:
        raise Exception('expected an even number of lines')
    samples = args.samples
    model_params = {
            "pa" : rat(25, 100),
            "pc" : rat(25, 100),
            "pg" : rat(25, 100),
            "pt" : rat(25, 100),
            "lambda" : rat(1, 1),
            "mu" : rat(2, 1),
            "tau" : rat(1, 10)}
    for i in range(npairs):
        print('benchmark', i, '...')
        a, b = lines[i*2 : (i+1)*2]
        d = dict(
            precision=args.precision,
            rtol=0.0,
            samples=samples,
            parameters=model_params,
            sequence_a=a,
            sequence_b=b)
        milliseconds = bench_pair(bench, d)
        print('seq1 length:', len(a), 'seq2 length:', len(b))
        print('elapsed milliseconds:', milliseconds / samples)
        print()


def pos_int(x):
    x = int(x)
    assert x > 0
    return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences', required=True,
            help='text file with raw sequence pairs')
    parser.add_argument('--samples', default=10, type=pos_int,
            help='number of runs per alignment for benchmarking')
    parser.add_argument('--precision', required=True,
            choices=('float', 'double', 'mag', 'arb256'))
    args = parser.parse_args()
    main(args)
