"""
Check whether or what proportion of alignments computed by arbtkf91-align
with float or double precision are optimal or canonical.

"""
from __future__ import print_function, division

from StringIO import StringIO
from subprocess import Popen, PIPE
import os
import json

from data_source import gen_files, gen_sequence_pairs


align = os.path.realpath('./bin/arbtkf91-align')
check = os.path.realpath('./bin/arbtkf91-check')


def align_pair(model_params, precision, a, b):
    d = dict(
        parameters=model_params,
        precision=precision,
        sequence_a=a,
        sequence_b=b)
    s_in = json.dumps(d)
    args = [align]
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    return json.loads(stdout_data)


def check_pair(model_params, a, b):
    d = dict(
        parameters=model_params,
        sequence_a=a,
        sequence_b=b)
    s_in = json.dumps(d)
    args = [check]
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    return json.loads(stdout_data)

def rat(a, b):
    return dict(num=a, denom=b)

def main():
    model_params = {
            "pa" : rat(27, 100),
            "pc" : rat(24, 100),
            "pg" : rat(26, 100),
            "pt" : rat(23, 100),
            "lambda" : rat(1, 1),
            "mu" : rat(2, 1),
            "tau" : rat(1, 10)}
    for name, fin in gen_files():
        print(name)
        sequence_pairs = list(gen_sequence_pairs(fin, force_acgt=True))
        for precision in 'float', 'double', 'mag':
            print('precision:', precision)
            ncanon = 0
            nopt = 0
            k = 0
            for a, b in sequence_pairs:
                d = align_pair(model_params, precision, a, b)
                a_aln = d['sequence_a']
                b_aln = d['sequence_b']
                js = check_pair(model_params, a_aln, b_aln)
                if js['alignment_is_canonical'] == 'yes':
                    ncanon += 1
                if js['alignment_is_optimal'] == 'yes':
                    nopt += 1
                k += 1
            print('total number of alignments:', k)
            print('number of optimal alignments:', nopt)
            print('number that are also canonical:', ncanon)
            print()


if __name__ == '__main__':
    main()
