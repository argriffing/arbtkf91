"""
Give arbtkf91-align some inputs that should cause failure.

"""
from __future__ import print_function, division

import os
import json
import copy
from subprocess import Popen, PIPE
from numpy.testing import assert_equal, assert_raises

align = os.path.realpath('./bin/arbtkf91-align')

root = {
        "precision" : "float",
        "rtol" : 0.0,
        "parameters" : {
            "pa" : {"num" : 25, "denom" : 100},
            "pc" : {"num" : 25, "denom" : 100},
            "pg" : {"num" : 25, "denom" : 100},
            "pt" : {"num" : 25, "denom" : 100},
            "lambda" : {"num" : 1, "denom" : 1},
            "mu" : {"num" : 2, "denom" : 1},
            "tau" : {"num" : 1, "denom" : 10}
        },
        "sequence_a" : "ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA",
        "sequence_b" : "AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA"}

class ReturnError(Exception):
    pass

def runjson(args, d):
    s_in = json.dumps(d)
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    out, err = data
    if p.returncode:
        raise ReturnError(err)
    else:
        return json.loads(out)

def test_ok():
    runjson([align], root)

def test_missing_precision():
    x = copy.deepcopy(root)
    del x['precision']
    assert_raises(ReturnError, runjson, [align], x)

def test_missing_parameters_pc():
    x = copy.deepcopy(root)
    del x['parameters']['pc']
    assert_raises(ReturnError, runjson, [align], x)

def test_missing_parameters_pc_denom():
    x = copy.deepcopy(root)
    del x['parameters']['pc']['denom']
    assert_raises(ReturnError, runjson, [align], x)
