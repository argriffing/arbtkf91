"""
Run programs that other people have written.

In this one, the third and fourth lines of output are the sequences.

"""
from __future__ import print_function, division

from StringIO import StringIO
from subprocess import Popen, PIPE
import os
import json

from data_source import gen_files, gen_sequence_pairs


align = os.path.realpath(os.path.join(
    '..', 'bioinf2015', 'implementation-team1',
    'build', 'mlalign'))

check = 'arbtkf91-check'


def align_pair(a, b):
    pairs = [
        ("runs", "1"),
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


def main():
    for name, fin in gen_files():
        ncanon = 0
        nopt = 0
        k = 0
        for a, b in gen_sequence_pairs(fin, force_acgt=True):
            out, err = align_pair(a, b)
            #print(out)
            #print(err)
            lines = [x.strip() for x in out.splitlines()]
            a_aln = lines[2]
            b_aln = lines[3]
            js = check_pair(a_aln, b_aln)
            if js['alignment_is_canonical'] == 'yes':
                ncanon += 1
            if js['alignment_is_optimal'] == 'yes':
                nopt += 1
            k += 1
        print(name)
        print('total number of alignments:', k)
        print('number of optimal alignments:', nopt)
        print('number that are also canonical:', ncanon)
        print()


if __name__ == '__main__':
    main()
