"""
Run programs that other people have written.

./tkf91 --sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA --sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA --lambda 1 --mu 2 --tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25

"""
from __future__ import print_function, division

import os
from subprocess import Popen, PIPE

from data_source import gen_files, gen_sequence_pairs


def align_pair(align, model_params, a, b):
    d = model_params[:]
    d.extend([
        ("sequence-1", a),
        ("sequence-2", b)])
    args = [align]
    for k, v in d:
        args.extend(['--'+k, v])
    print(args)
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    return p.communicate()


def main():
    align = os.path.realpath('../../stamatakis/tkf91_scaling/tkf91')
    model_params = [
            ("pa", "0.27"),
            ("pc", "0.24"),
            ("pg", "0.26"),
            ("pt", "0.23"),
            ("lambda", "1"),
            ("mu", "2"),
            ("tau", "0.1"),
            ]
    #a = 'ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA'
    #b = 'AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA'
    #a = a*10
    #b = b*10
    # 45 ok
    # 46 fail
    # 47
    # 48
    # 49
    # 50 fail
    n = 45
    a = 'A'*n
    b = 'C'*n
    out, err = align_pair(align, model_params, a, b)
    if 'failed' in err:
        print('it seemed to fail with the error message')
        print(err)
    else:
        print('it seemed to work')
    return
    for name, fin in gen_files():
        print(name)
        with open("out-" + name, "wt") as fout:
            with open("err-" + name, "wt") as ferr:
                for a, b in gen_sequence_pairs(fin, force_acgt=True):
                    out, err = align_pair(align, model_params, a, b)
                    print(out, file=fout)
                    print(err, file=ferr)
                    fout.flush()
                    ferr.flush()


if __name__ == '__main__':
    main()
