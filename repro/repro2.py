"""
Run programs that other people have written.

For this one, it looks like the output can be parsed by looking for the line
that says 'log_caching_round_up' and then the next two lines are the
nucleotide sequences.

"""
from __future__ import print_function, division

from subprocess import Popen, PIPE
import os
import json
import argparse

from data_source import gen_files, gen_sequence_pairs


check = 'arbtkf91-check'


def align_pair(align, a, b):
    pairs = [
        ("output-alignment", "1"),
        ("pa", "0.27"),
        ("pc", "0.24"),
        ("pg", "0.26"),
        ("pt", "0.23"),
        ("lambda", "1"),
        ("mu", "2"),
        ("tau", "0.1"),
        ("sequence-1", a),
        ("sequence-2", b),
        ]
    args = [align]
    for k, v in pairs:
        args.extend(['--'+k, v])
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    return p.communicate()


def runjson(args, d):
    s_in = json.dumps(d)
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    return json.loads(stdout_data)

def rat(a, b):
    return dict(num=a, denom=b)

def check_pair(a, b):
    model_params = {
            "pa" : rat(27, 100),
            "pc" : rat(24, 100),
            "pg" : rat(26, 100),
            "pt" : rat(23, 100),
            "lambda" : rat(1, 1),
            "mu" : rat(2, 1),
            "tau" : rat(1, 10)}
    j_in = dict(
        parameters=model_params,
        sequence_a=a,
        sequence_b=b)
    return runjson([check], j_in)


def main(args):
    align = os.path.realpath(args.exe)
    data_path = os.path.realpath(args.bench_data)
    for name, fin in gen_files(data_path):
        ncanon = 0
        nopt = 0
        k = 0
        for a, b in gen_sequence_pairs(fin, force_acgt=True):
            out, err = align_pair(align, a, b)
            lines = [x.strip() for x in out.splitlines()]
            nlines = len(lines)
            found = False
            for i in range(nlines):
                if lines[i] == 'log_caching_round_up':
                    a_aln = lines[i+1]
                    b_aln = lines[i+2]
                    found = True
                    break
            if not found:
                raise Exception('failed to understand the output')
            js = check_pair(a_aln, b_aln)
            if js['alignment_is_canonical']:
                ncanon += 1
            if js['alignment_is_optimal']:
                nopt += 1
            k += 1
        print(name)
        print('total number of alignments:', k)
        print('number of optimal alignments:', nopt)
        print('number that are also canonical:', ncanon)
        print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bench-data', required=True,
            help='benchmark data directory')
    parser.add_argument('--exe', required=True,
            help='TKFLOG_CACHING_ROUND_UP alignment executable')
    args = parser.parse_args()
    main(args)
