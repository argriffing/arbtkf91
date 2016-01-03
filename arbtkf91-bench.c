/*
 * Benchmark an alignment strategy.
 * Input and output uses json.
 *
 * input:
 * {
 * "precision" : "float" | "double" | "exact",
 * "samples" : integer,
 * "pa_n" : integer, "pa_d" : integer,
 * "pc_n" : integer, "pc_d" : integer,
 * "pg_n" : integer, "pg_d" : integer,
 * "pt_n" : integer, "pt_d" : integer,
 * "lambda_n" : integer, "lambda_d" : integer,
 * "mu_n" : integer, "mu_d" : integer,
 * "tau_n" : integer, "tau_d" : integer,
 * "sequence_a" : string,
 * "sequence_b" : string
 * }
 *
 * output:
 * {
 * "ticks_per_second" : integer,
 * "elapsed_ticks" : [integer, integer, ..., integer],
 * "sequence_a" : string,
 * "sequence_b" : string,
 * "verified" : bool
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
#include "tkf91_dp_bound.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"
#include "model_params.h"



void
solve(tkf91_dp_fn f, solution_t sol, const model_params_t p,
        const slong *A, slong len_A, const slong *B, slong len_B);


json_t *run(void * userdata, json_t *j_in);

json_t *run(void * userdata, json_t *j_in)
{
    json_t *args;
    json_t *j_out;
    model_params_t p;
    slong len_A, len_B;
    slong *A;
    slong *B;
    solution_t sol;

    model_params_init(p);

    if (userdata)
    {
        fprintf(stderr, "error: unexpected userdata\n");
        abort();
    }

    args = j_in;

    /* read the model parameter values */
    _json_object_get_fmpq(p->pi+0, args, "pa_n", "pa_d");
    _json_object_get_fmpq(p->pi+1, args, "pc_n", "pc_d");
    _json_object_get_fmpq(p->pi+2, args, "pg_n", "pg_d");
    _json_object_get_fmpq(p->pi+3, args, "pt_n", "pt_d");
    _json_object_get_fmpq(p->lambda, args, "lambda_n", "lambda_d");
    _json_object_get_fmpq(p->mu, args, "mu_n", "mu_d");
    _json_object_get_fmpq(p->tau, args, "tau_n", "tau_d");

    /* read the two unaligned sequences */
    A = _json_object_get_sequence(&len_A, args, "sequence_a");
    B = _json_object_get_sequence(&len_B, args, "sequence_b");

    /* read the requested precision and the requested number of samples */
    const char * precision_str;
    slong samples;
    precision_str = _json_object_get_string(args, "precision");
    samples = _json_object_get_si(args, "samples");
    if (samples < 1)
    {
        printf("at least one sample is required\n");
        abort();
    }

    /* dispatch */
    tkf91_dp_fn f = NULL;
    if (strcmp(precision_str, "float") == 0)
    {
        f = tkf91_dp_f;
    }
    else if (strcmp(precision_str, "double") == 0)
    {
        f = tkf91_dp_d;
    }
    else if (strcmp(precision_str, "exact") == 0)
    {
        ;
        f = tkf91_dp_bound;
    }
    else
    {
        printf("expected the precision string to be one of ");
        printf("{'float' | 'double' | 'exact'}\n");
        abort();
    }

    int i;
    clock_t start, diff;
    json_t *elapsed_ticks;
    elapsed_ticks = json_array();
    for (i = 0; i < samples; i++)
    {
        start = clock();
        solution_init(sol, len_A + len_B);
        solve(f, sol, p, A, len_A, B, len_B);
        if (i < samples - 1)
        {
            solution_clear(sol);
        }
        diff = clock() - start;
        json_array_append_new(elapsed_ticks, json_integer((json_int_t) diff));
    }

    j_out = json_pack("{s:i, s:o, s:s, s:s, s:b}",
            "ticks_per_second", (json_int_t) CLOCKS_PER_SEC,
            "elapsed_ticks", elapsed_ticks,
            "sequence_a", sol->A,
            "sequence_b", sol->B,
            "verified", sol->optimality_flag);

    flint_free(A);
    flint_free(B);
    solution_clear(sol);
    model_params_clear(p);

    return j_out;
}



void
solve(tkf91_dp_fn f, solution_t sol, const model_params_t p,
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
    req->png_filename = NULL;
    req->trace = 1;

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
