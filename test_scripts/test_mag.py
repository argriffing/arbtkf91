"""
Run arbtkf91-align on random inputs and check the output with arbtkf91-check.

"""
from __future__ import print_function, division

from subprocess import Popen, PIPE
import random
import os
import json
from numpy.testing import assert_equal

align = 'arbtkf91-align'
check = 'arbtkf91-check'

def runjson(args, d):
    s_in = json.dumps(d)
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    data = p.communicate(input=s_in)
    out, err = data
    try:
        return json.loads(out)
    except ValueError as e:
        print(out)
        print(err)
        raise

def rat(a, b):
    return dict(num=a, denom=b)


def sample_params():

    # Sample the nucleotide distribution.
    # Require that the probabilities are positive and sum to 1.
    pi_num = [random.randrange(1, 100) for i in range(4)]
    pi_denom = sum(pi_num)

    # Sample parameters of the birth-death process.
    # Require that 0 < lambda < mu.
    lambda_num = random.randrange(1, 100)
    mu_num = lambda_num + random.randrange(1, 100)
    bd_denom = random.randrange(1, 100)

    model_params = {
            "pa" : rat(pi_num[0], pi_denom),
            "pc" : rat(pi_num[1], pi_denom),
            "pg" : rat(pi_num[2], pi_denom),
            "pt" : rat(pi_num[3], pi_denom),
            "lambda" : rat(lambda_num, bd_denom),
            "mu" : rat(mu_num, bd_denom),
            "tau" : rat(random.randrange(1, 10), random.randrange(1, 100))}

    return model_params


def sample_sequence():
    ACGT = 'ACGT'
    alpha = random.choice((5, 50, 500))
    n = random.randrange(2, alpha)
    return ''.join(random.choice(ACGT) for i in range(n))

def sample_sequences():
    a = sample_sequence()
    b = sample_sequence()
    return a, b

def check_mag(model_params, a, b):
    j_in = dict(
        parameters=model_params,
        rtol=0.0,
        precision='mag',
        sequence_a=a,
        sequence_b=b)
    d = runjson([align], j_in)
    j_in = dict(
        parameters=model_params,
        sequence_a=d['sequence_a'],
        sequence_b=d['sequence_b'])
    d = runjson([check], j_in)
    assert_equal(d['alignment_is_canonical'], True)
    assert_equal(d['alignment_is_optimal'], True)

def check_for_smoke(precision, rtol, model_params, a, b):
    # smoke test for float and double precision
    j_in = dict(
        parameters=model_params,
        rtol=rtol,
        precision=precision,
        sequence_a=a,
        sequence_b=b)
    d = runjson([align], j_in)
    j_in = dict(
        parameters=model_params,
        sequence_a=d['sequence_a'],
        sequence_b=d['sequence_b'])
    d = runjson([check], j_in)

def test_mag():
    random.seed(1234)
    nsamples = 20
    for i in range(nsamples):
        model_params = sample_params()
        a, b = sample_sequences()
        check_mag(model_params, a, b)

def test_smoke_float():
    random.seed(1234)
    nsamples = 20
    precision = 'float'
    for rtol in 0.0, 1e-2, 1e-7:
        for i in range(nsamples):
            model_params = sample_params()
            a, b = sample_sequences()
            check_for_smoke(precision, rtol, model_params, a, b)

def test_smoke_double():
    random.seed(1234)
    nsamples = 20
    precision = 'double'
    for rtol in 0.0, 1e-2, 1e-7:
        for i in range(nsamples):
            model_params = sample_params()
            a, b = sample_sequences()
            check_for_smoke(precision, rtol, model_params, a, b)
