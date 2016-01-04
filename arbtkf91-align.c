/*
 * Align sequences.
 * Input and output uses json.
 * The output has all of the terms expected by the 'arbtkf91-check' tool,
 * and it also has a "verified" boolean.
 *
 * input:
 * {
 * "precision" : "float" | "double" | "high", "exact",
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
 * "pa_n" : integer, "pa_d" : integer,
 * "pc_n" : integer, "pc_d" : integer,
 * "pg_n" : integer, "pg_d" : integer,
 * "pt_n" : integer, "pt_d" : integer,
 * "lambda_n" : integer, "lambda_d" : integer,
 * "mu_n" : integer, "mu_d" : integer,
 * "tau_n" : integer, "tau_d" : integer,
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
#include "tkf91_dp_r.h"
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

    /* read the requested precision */
    const char * precision_str;
    precision_str = _json_object_get_string(args, "precision");

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
    else if (strcmp(precision_str, "high") == 0)
    {
        f = tkf91_dp_r;
    }
    else if (strcmp(precision_str, "exact") == 0)
    {
        f = tkf91_dp_bound;
    }
    else
    {
        printf("expected the precision string to be one of ");
        printf("{\"float\" | \"double\" | \"exact\" | \"high\"}\n");
        abort();
    }

    solution_init(sol, len_A + len_B);
    solve(f, sol, p, A, len_A, B, len_B);

    /*
     * Copy the input json object and delete a few keys.
     * The remaining keys are kept so that the output json object
     * can be passed directly to 'arbtkf91-check'.
     */
    j_out = json_copy(j_in);
    json_object_del(j_out, "precision");
    json_object_del(j_out, "sequence_a");
    json_object_del(j_out, "sequence_b");

    /* put the aligned sequences into the output object */
    json_object_set_new(j_out, "sequence_a", json_string(sol->A));
    json_object_set_new(j_out, "sequence_b", json_string(sol->B));
    json_object_set_new(j_out, "verified", json_boolean(sol->optimality_flag));

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