"""
Check what happens when bad arguments are passed the the programs.

"""
from __future__ import print_function, division

import os
import json
import copy
from subprocess import Popen, PIPE
from numpy.testing import assert_equal, assert_raises, TestCase

align = 'arbtkf91-align'
check = 'arbtkf91-check'
bench = 'arbtkf91-bench'
count = 'arbtkf91-count'

too_small_pi = {
    "pa" : {"num" : 24, "denom" : 100},
    "pc" : {"num" : 25, "denom" : 100},
    "pg" : {"num" : 25, "denom" : 100},
    "pt" : {"num" : 25, "denom" : 100},
    }

too_large_pi = {
    "pa" : {"num" : 26, "denom" : 100},
    "pc" : {"num" : 25, "denom" : 100},
    "pg" : {"num" : 25, "denom" : 100},
    "pt" : {"num" : 25, "denom" : 100},
    }

bench_root = {
        "precision" : "float",
        "rtol" : 0.0,
        "samples" : 10,
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

count_root = {
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


class CheckArgs(object):
    args = None
    root = None

    def setUp(self):
        self.x = copy.deepcopy(self.root)

    def check_return_error(self):
        assert_raises(ReturnError, runjson, self.args, self.x)

    def test_ok(self):
        runjson(self.args, self.root)

    def test_too_small_pi_sum(self):
        self.x['parameters'].update(too_small_pi)
        self.check_return_error()

    def test_too_large_pi_sum(self):
        self.x['parameters'].update(too_large_pi)
        self.check_return_error()

    def test_negative_numerator(self):
        for p in 'pa', 'pc', 'pg', 'pt', 'lambda', 'mu', 'tau':
            self.x = copy.deepcopy(self.root)
            self.x['parameters'][p]['num'] *= -1
            self.check_return_error()

    def test_negative_denominator(self):
        for p in 'pa', 'pc', 'pg', 'pt', 'lambda', 'mu', 'tau':
            self.x = copy.deepcopy(self.root)
            self.x['parameters'][p]['denom'] *= -1
            self.check_return_error()

    def test_zero_numerator(self):
        for p in 'pa', 'pc', 'pg', 'pt', 'lambda', 'mu', 'tau':
            self.x = copy.deepcopy(self.root)
            self.x['parameters'][p]['num'] = 0
            self.check_return_error()

    def test_zero_denominator(self):
        for p in 'pa', 'pc', 'pg', 'pt', 'lambda', 'mu', 'tau':
            self.x = copy.deepcopy(self.root)
            self.x['parameters'][p]['denom'] = 0
            self.check_return_error()

    def test_mu_equals_lambda(self):
        mu = self.x['parameters']['mu']
        self.x['parameters']['mu'] = copy.deepcopy(mu)
        self.x['parameters']['lambda'] = copy.deepcopy(mu)
        self.check_return_error()

    def test_mu_less_than_lambda(self):
        mu = self.x['parameters']['mu']
        self.x['parameters']['mu'] = copy.deepcopy(mu)
        self.x['parameters']['mu']['denom'] *= 2
        self.x['parameters']['lambda'] = copy.deepcopy(mu)
        self.check_return_error()

    def test_missing_parameters(self):
        del self.x['parameters']
        self.check_return_error()

    def test_missing_parameters_pc(self):
        del self.x['parameters']['pc']
        self.check_return_error()

    def test_missing_parameters_pc_denom(self):
        del self.x['parameters']['pc']['denom']
        self.check_return_error()

    def test_missing_sequence(self):
        del self.x['sequence_a']
        self.check_return_error()



class TestAlignArgs(CheckArgs, TestCase):
    root = align_root
    args = [align]


class TestCheckArgs(CheckArgs, TestCase):
    root = check_root
    args = [check]


class TestBenchArgs(CheckArgs, TestCase):
    root = bench_root
    args = [bench]

class TestCountArgs(CheckArgs, TestCase):
    root = count_root
    args = [count]
