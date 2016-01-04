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

check = os.path.realpath('./bin/arbtkf91-check')


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


def check_pair(a, b):
    model_params = dict(
            pa_n=27, pa_d=100,
            pc_n=24, pc_d=100,
            pg_n=26, pg_d=100,
            pt_n=23, pt_d=100,
            lambda_n=1, lambda_d=1,
            mu_n=2, mu_d=1,
            tau_n=1, tau_d=10)
    d = model_params.copy()
    d.update(sequence_a=a, sequence_b=b)
    s_in = json.dumps(d)
    args = [check]
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    return json.loads(stdout_data)


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
