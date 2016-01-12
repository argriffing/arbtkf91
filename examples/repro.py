"""
The reference alignment program does not work for me...

"""
from __future__ import print_function, division

import argparse
import os
from subprocess import Popen, PIPE


def align_pair(align, model_params, a, b):
    d = model_params[:]
    d.extend([
        ("sequence-1", a),
        ("sequence-2", b)])
    args = [align]
    for k, v in d:
        args.extend(['--'+k, v])
    print('/path/to/exe', ' '.join(args[1:]))
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    return p.communicate()


def run(align, model_params, a, b):
    out, err = align_pair(align, model_params, a, b)
    if 'failed' in err:
        print('It seemed to fail with the following error message:')
        print(err)
    else:
        print('It seemed to work.')


def main(args):
    align = os.path.realpath(args.ref_exe)
    model_params = [
            ("pa", "0.27"),
            ("pc", "0.24"),
            ("pg", "0.26"),
            ("pt", "0.23"),
            ("lambda", "1"),
            ("mu", "2"),
            ("tau", "0.1"),
            ]

    print('first example...')
    a = 'ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA'
    b = 'AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA'
    run(align, model_params, a, b)
    print()

    print('second example...')
    n = 45
    a = 'A'*n
    b = 'C'*n
    run(align, model_params, a, b)
    print()

    print('third example...')
    n = 46
    a = 'A'*n
    b = 'C'*n
    run(align, model_params, a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-exe', required=True,
            help='reference alignment executable')
    args = parser.parse_args()
    main(args)
