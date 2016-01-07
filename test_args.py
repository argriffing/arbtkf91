"""
Check what happens when bad arguments are passed
to programs like arbtkf91-align and arbtkf91-check.

"""
from __future__ import print_function, division

import os
import json
import copy
from subprocess import Popen, PIPE
from numpy.testing import assert_equal, assert_raises, TestCase

align = os.path.realpath('./bin/arbtkf91-align')
check = os.path.realpath('./bin/arbtkf91-check')

align_root = {
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

check_root = {
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
        "sequence_b" : "AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA---"}

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


class TestAlignArgs(TestCase):

    def test_ok(self):
        runjson([align], align_root)

    def test_missing_precision(self):
        x = copy.deepcopy(align_root)
        del x['precision']
        assert_raises(ReturnError, runjson, [align], x)

    def test_missing_parameters_pc(self):
        x = copy.deepcopy(align_root)
        del x['parameters']['pc']
        assert_raises(ReturnError, runjson, [align], x)

    def test_missing_parameters_pc_denom(self):
        x = copy.deepcopy(align_root)
        del x['parameters']['pc']['denom']
        assert_raises(ReturnError, runjson, [align], x)


class TestCheckArgs(TestCase):

    def test_ok(self):
        runjson([check], check_root)

    def test_missing_sequence(self):
        x = copy.deepcopy(check_root)
        del x['sequence_a']
        assert_raises(ReturnError, runjson, [check], x)

    def test_missing_parameters_pc(self):
        x = copy.deepcopy(check_root)
        del x['parameters']['pc']
        assert_raises(ReturnError, runjson, [check], x)

    def test_missing_parameters_pc_denom(self):
        x = copy.deepcopy(check_root)
        del x['parameters']['pc']['denom']
        assert_raises(ReturnError, runjson, [check], x)
