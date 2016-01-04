"""
yield stuff from a hardcoded data source

"""
from __future__ import print_function, division

import os

import numpy as np
from numpy.testing import assert_array_equal, assert_equal
from Bio import SeqIO

__all__ = ['gen_files', 'gen_sequence_pairs']

mypath = os.path.realpath('../../stamatakis/benchMark_data')

def _normalized_seq(s):
    return ''.join(_normalized_chr(c) for c in s)

def _normalized_chr(c):
    if c in 'ACGT':
        return c
    elif c.isalpha:
        return 'A'
    else:
        msg = ('weird character:', c)
        raise Exception(msg)


def gen_sequence_pairs(fin, force_acgt=False):
    # yield (10 choose 2) = 45 nucleotide sequence pairs
    fasta_objects = list(SeqIO.parse(fin, 'fasta'))
    sequences = [str(x.seq) for x in fasta_objects]
    available = len(sequences)
    requested = 10
    indices = _select_indices(available, requested)
    selection = [sequences[i] for i in indices]
    assert_equal(len(selection), requested)
    k = 0
    for i in range(requested):
        for j in range(i):
            a = selection[i]
            b = selection[j]
            if force_acgt:
                a = _normalized_seq(a)
                b = _normalized_seq(b)
            yield a, b
            k += 1
    assert_equal(k, 45)

def gen_files():
    # yield (name, handle) pairs
    for filename in os.listdir(mypath):
        if 'unaligned' in filename:
            fullpath = os.path.join(mypath, filename)
            with open(fullpath) as fin:
                yield filename, fin


def _select_indices(available, requested):
    incr = available // requested
    return [i*incr for i in range(requested)]


def test():
    indices = _select_indices(60, 10)
    assert_array_equal(indices[:3], [0, 6, 12])
    assert_array_equal(indices[-1:], [54])


if __name__ == '__main__':
    test()
