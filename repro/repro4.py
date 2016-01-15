"""
Check whether or what proportion of alignments computed by arbtkf91-align
with float or double precision are optimal or canonical.

"""
from __future__ import print_function, division

from subprocess import Popen, PIPE
import argparse
import os
import json

from data_source import gen_files, gen_sequence_pairs

align = 'arbtkf91-align'
check = 'arbtkf91-check'

def runjson(args, d):
    s_in = json.dumps(d)
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    return json.loads(stdout_data)

def rat(a, b):
    return dict(num=a, denom=b)


def main(args):
    data_path = os.path.realpath(args.bench_data)
    model_params = {
            "pa" : rat(27, 100),
            "pc" : rat(24, 100),
            "pg" : rat(26, 100),
            "pt" : rat(23, 100),
            "lambda" : rat(1, 1),
            "mu" : rat(2, 1),
            "tau" : rat(1, 10)}
    for name, fin in gen_files(data_path):
        print(name)
        print('precision={} rtol={}'.format(args.precision, args.rtol))
        sequence_pairs = list(gen_sequence_pairs(fin, force_acgt=True))
        ncanon = 0
        nopt = 0
        k = 0
        for a, b in sequence_pairs:
            j_in = dict(
                parameters=model_params,
                rtol=args.rtol,
                precision=args.precision,
                sequence_a=a,
                sequence_b=b)
            d = runjson([align], j_in)
            j_in = dict(
                parameters=model_params,
                sequence_a=d['sequence_a'],
                sequence_b=d['sequence_b'])
            d = runjson([check], j_in)
            if d['alignment_is_canonical']:
                ncanon += 1
            if d['alignment_is_optimal']:
                nopt += 1
            k += 1
        print('total number of alignments:', k)
        print('number of optimal alignments:', nopt)
        print('number that are also canonical:', ncanon)
        print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bench-data', required=True,
            help='benchmark data directory')
    parser.add_argument('--precision', required=True,
            choices=('float', 'double', 'mag', 'high'))
    parser.add_argument('--rtol', type=float, default=0.0,
            help="relative tolerance for float and double precision")
    args = parser.parse_args()
    main(args)
