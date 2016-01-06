from __future__ import print_function, division

import os
from subprocess import Popen, PIPE
import json

import data_source

def rat(a, b):
    return dict(num=a, denom=b)

def bench_pair(bench, model_params, a, b):
    d = dict(
        precision='float',
        samples=10,
        parameters=model_params,
        sequence_a=a,
        sequence_b=b)
    s_in = json.dumps(d)
    args = [bench]
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    stdout_data, stderr_data = data
    out = json.loads(stdout_data)
    elapsed = out['elapsed_ticks']
    print(elapsed)
    return sum(elapsed)


def main():
    bench = os.path.realpath('./bin/arbtkf91-bench')
    model_params = {
            "pa" : rat(27, 100),
            "pc" : rat(24, 100),
            "pg" : rat(26, 100),
            "pt" : rat(23, 100),
            "lambda" : rat(1, 1),
            "mu" : rat(2, 1),
            "tau" : rat(1, 10)}
    for name, fin in data_source.gen_files():
        print(name)
        total_ticks = 0
        for a, b in data_source.gen_sequence_pairs(fin):
            total_ticks += bench_pair(bench, model_params, a, b)
        print(total_ticks / 10)
        print()


if __name__ == '__main__':
    main()
