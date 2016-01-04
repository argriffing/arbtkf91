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


def gen_sequence_pairs(fin):
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
