/*
 * Align sequences.
 * Input and output uses json.
 * The output has all of the terms expected by the 'arbtkf91-check' tool,
 * and it also has a "verified" boolean.
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
    json_t *parameters;
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

    flags = JSON_STRICT;
    result = json_unpack_ex(root, &err, flags,
            "{s:O, s:s, s:s, s:F, s:s}",
            "parameters", &parameters,
            "sequence_a", &sequence_a,
            "sequence_b", &sequence_b,
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
    solution_init(sol, len_A + len_B);

    /* dispatch, initializing a tableau if necessary */
    dp_mat_t tableau;
    tkf91_dp_fn f = NULL;
    if (strcmp(precision, "float") == 0) {
        f = tkf91_dp_f;
    }
    else if (strcmp(precision, "double") == 0) {
        f = tkf91_dp_d;
    }
    else if (strcmp(precision, "mag") == 0) {
        f = tkf91_dp_bound;
        dp_mat_init(tableau, nrows, ncols);
        sol->mat = tableau;
    }
    else if (strcmp(precision, "arb256") == 0) {
        f = tkf91_dp_r;
        dp_mat_init(tableau, nrows, ncols);
        sol->mat = tableau;
    }
    else {
        printf("expected the precision string to be one of ");
        printf("{float | double | mag | arb256}\n");
        abort();
    }

    solve(f, sol, rtol, p, A, len_A, B, len_B);

    j_out = json_pack("{s:o, s:s, s:s, s:b}",
            "parameters", parameters,
            "sequence_a", sol->A,
            "sequence_b", sol->B,
            "verified", sol->optimality_flag);

    flint_free(A);
    flint_free(B);
    solution_clear(sol);
    model_params_clear(p);
    if (sol->mat)
    {
        dp_mat_clear(sol->mat);
    }

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
    req->png_filename = NULL;
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
