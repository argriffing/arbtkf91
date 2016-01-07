/*
 * "pa_n" : integer, "pa_d" : integer,
 * "pc_n" : integer, "pc_d" : integer,
 * "pg_n" : integer, "pg_d" : integer,
 * "pt_n" : integer, "pt_d" : integer,
 * "lambda_n" : integer, "lambda_d" : integer,
 * "mu_n" : integer, "mu_d" : integer,
 * "tau_n" : integer, "tau_d" : integer,
 *
 * "parameters" : {
 * "pa" : {"num" : integer, "denom" : integer},
 * "pc" : {"num" : integer, "denom" : integer},
 * "pg" : {"num" : integer, "denom" : integer},
 * "pt" : {"num" : integer, "denom" : integer},
 * "lambda" : {"num" : integer, "denom" : integer},
 * "mu" : {"num" : integer, "denom" : integer},
 * "tau" : {"num" : integer, "denom" : integer}
 * }
 */

#include "flint/fmpq.h"

#include "jsonutil.h"
#include "model_params.h"
#include "json_model_params.h"

int
_json_get_fmpq(fmpq_t x, json_t * root)
{
    json_error_t error;
    size_t flags;
    flags = 0;
    return _json_get_fmpq_ex(x, root, &error, flags);
}

int
_json_get_fmpq_ex(fmpq_t x, json_t * root, json_error_t *perror, size_t flags)
{
    int rnum, rdenom, result;
    slong num;
    ulong denom;
    result = json_unpack_ex(root, perror, flags, "{s:i, s:i}",
            "num", &rnum,
            "denom", &rdenom);
    if (result)
    {
        return result;
    }
    if (rdenom < 0)
    {
        rnum = -rnum;
        rdenom = -rdenom;
    }
    num = (slong) rnum;
    denom = (ulong) rdenom;
    fmpq_set_si(x, num, denom);
    return 0;
}

int
_json_get_model_params(model_params_t p, json_t * root)
{
    json_error_t error;
    size_t flags;
    flags = 0;
    return _json_get_model_params_ex(p, root, &error, flags);
}

int
_json_get_model_params_ex(model_params_t p, json_t * root,
        json_error_t *e, size_t flags)
{
    int r;
    json_t *pa, *pc, *pg, *pt, *lambda, *mu, *tau;

    r = json_unpack_ex(root, e, flags,
            "{s:o, s:o, s:o, s:o, s:o, s:o, s:o}",
            "pa", &pa,
            "pc", &pc,
            "pg", &pg,
            "pt", &pt,
            "lambda", &lambda,
            "mu", &mu,
            "tau", &tau);
    if (r) return r;

    r = _json_get_fmpq_ex(p->pi + 0, pa, e, flags); if (r) return r;
    r = _json_get_fmpq_ex(p->pi + 1, pc, e, flags); if (r) return r;
    r = _json_get_fmpq_ex(p->pi + 2, pg, e, flags); if (r) return r;
    r = _json_get_fmpq_ex(p->pi + 3, pt, e, flags); if (r) return r;
    r = _json_get_fmpq_ex(p->lambda, lambda, e, flags); if (r) return r;
    r = _json_get_fmpq_ex(p->mu, mu, e, flags); if (r) return r;
    r = _json_get_fmpq_ex(p->tau, tau, e, flags); if (r) return r;
    return 0;
}
