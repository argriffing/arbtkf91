/*
 * Benchmark an alignment strategy.
 * The json format is used for both the input and the output.
 *
 * output:
 * {
 * "ticks_per_second" : integer,
 * "elapsed_ticks" : [integer, integer, ..., integer],
 * "sequence_a" : string,
 * "sequence_b" : string
 * }
 *
 */

#include <time.h>

#include "flint/flint.h"
#include "flint/fmpq.h"

#include "jansson.h"

#include "runjson.h"
#include "jsonutil.h"
#include "tkf91_dp_f.h"
#include "tkf91_dp_d.h"
#include "tkf91_dp_r.h"
#include "tkf91_dp_bound.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"
#include "model_params.h"
#include "json_model_params.h"



void
solve(tkf91_dp_fn f, solution_t sol, double rtol, const model_params_t p,
        const slong *A, slong len_A, const slong *B, slong len_B);


json_t *run(void * userdata, json_t *root);

json_t *run(void * userdata, json_t *root)
{
    json_t *j_out;
    model_params_t p;
    slong len_A, len_B;
    slong *A;
    slong *B;
    solution_t sol;
    int result;
    int samples;
    json_t * parameters;
    const char * sequence_a;
    const char * sequence_b;
    const char * precision;
    double rtol;
    json_error_t err;
    size_t flags;
    slong nrows, ncols;

    if (userdata)
    {
        fprintf(stderr, "error: unexpected userdata\n");
        abort();
    }

    /* default values of optional json arguments */
    rtol = 0;
    precision = NULL;

    flags = JSON_STRICT;
    result = json_unpack_ex(root, &err, flags,
            "{s:o, s:s, s:s, s:i, s?F, s?s}",
            "parameters", &parameters,
            "sequence_a", &sequence_a,
            "sequence_b", &sequence_b,
            "samples", &samples,
            "rtol", &rtol,
            "precision", &precision);
    if (result)
    {
        fprintf(stderr, "error: on line %d: %s\n", err.line, err.text);
        abort();
    }

    model_params_init(p);
    result = _json_get_model_params_ex(p, parameters, &err, flags);
    if (result)
    {
        fprintf(stderr, "error: on line %d: %s\n", err.line, err.text);
        abort();
    }
    result = model_params_validate(p);
    if (result)
    {
        fprintf(stderr, "invalid model parameters\n");
        abort();
    }

    /* read the two unaligned sequences */

    len_A = strlen(sequence_a);
    A = flint_malloc(len_A * sizeof(slong));
    _fill_sequence_vector(A, sequence_a, len_A);

    len_B = strlen(sequence_b);
    B = flint_malloc(len_B * sizeof(slong));
    _fill_sequence_vector(B, sequence_b, len_B);

    nrows = len_A + 1;
    ncols = len_B + 1;

    /* dispatch */
    tkf91_dp_fn f = NULL;
    int requires_tableau = 0;
    if (precision == NULL) {
        f = tkf91_dp_high;
        requires_tableau = 1;
    }
    else if (strcmp(precision, "float") == 0) {
        f = tkf91_dp_f;
    }
    else if (strcmp(precision, "double") == 0) {
        f = tkf91_dp_d;
    }
    else if (strcmp(precision, "mag") == 0) {
        f = tkf91_dp_mag;
        requires_tableau = 1;
    }
    else if (strcmp(precision, "high") == 0) {
        f = tkf91_dp_high;
        requires_tableau = 1;
    }
    else
    {
        fprintf(stderr, "expected the precision string to be one of ");
        fprintf(stderr, "{'float' | 'double' | 'mag' | 'high'}\n");
        abort();
    }

    int i;
    clock_t start, diff;
    json_t *elapsed_ticks;
    elapsed_ticks = json_array();
    dp_mat_t tableau;
    for (i = 0; i < samples; i++)
    {
        start = clock();
        solution_init(sol, len_A + len_B);
        if (requires_tableau)
        {
            dp_mat_init(tableau, nrows, ncols);
            sol->mat = tableau;
        }
        solve(f, sol, rtol, p, A, len_A, B, len_B);
        if (requires_tableau)
        {
            dp_mat_clear(tableau);
        }
        if (i < samples - 1)
        {
            solution_clear(sol);
        }
        diff = clock() - start;
        json_array_append_new(elapsed_ticks, json_integer((json_int_t) diff));
    }

    j_out = json_pack("{s:i, s:o, s:s, s:s}",
            "ticks_per_second", (json_int_t) CLOCKS_PER_SEC,
            "elapsed_ticks", elapsed_ticks,
            "sequence_a", sol->A,
            "sequence_b", sol->B);

    flint_free(A);
    flint_free(B);
    solution_clear(sol);
    model_params_clear(p);

    return j_out;
}



void
solve(tkf91_dp_fn f, solution_t sol, double rtol, const model_params_t p,
        const slong *A, slong szA, const slong *B, slong szB)
{
    tkf91_rationals_t r;
    tkf91_expressions_t expressions;
    tkf91_generator_indices_t generators;
    fmpz_mat_t mat;
    expr_ptr * expressions_table;
    request_t req;

    /* expressions registry and (refining) generator registry */
    reg_t er;
    rgen_reg_ptr gr;

    reg_init(er);
    tkf91_rationals_init(r, p->lambda, p->mu, p->tau, p->pi);
    tkf91_expressions_init(expressions, er, r);

    gr = rgen_reg_new();
    tkf91_rgenerators_init(generators, gr, r, expressions, A, szA, B, szB);
    rgen_reg_finalize(gr, er);
    fmpz_mat_init(mat, rgen_reg_nrows(gr), rgen_reg_ncols(gr));
    rgen_reg_get_matrix(mat, gr);

    rgen_reg_clear(gr);
    tkf91_rationals_clear(r);

    expressions_table = reg_vec(er);

    /* init request object */
    req->trace = 1;
    req->rtol = rtol;

    f(sol, req, mat, expressions_table, generators, A, szA, B, szB);

    fmpz_mat_clear(mat);
    flint_free(expressions_table);

    reg_clear(er);
    tkf91_expressions_clear(expressions);
}



int main(void)
{
    json_hom_t hom;
    hom->userdata = NULL;
    hom->clear = NULL;
    hom->f = run;
    int result = run_json_script(hom);

    flint_cleanup();
    return result;
}
