/*
 * Count the number of optimal sequences.
 * Input and output uses json for consistency with the other scripts,
 * although that seems somewhat silly in this case because the
 * output is just a single number that is likely so large that it
 * cannot even be represented in json except as a string.
 */

#include "flint/flint.h"
#include "flint/fmpq.h"

#include "jansson.h"

#include "runjson.h"
#include "jsonutil.h"
#include "tkf91_dp_r.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"
#include "model_params.h"
#include "json_model_params.h"
#include "count_solutions.h"


void solve(fmpz_t res, solution_t sol, const model_params_t p,
        const slong *A, slong szA, const slong *B, slong szB);


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
    json_error_t err;
    size_t flags;
    slong nrows, ncols;
    dp_mat_t tableau;

    if (userdata)
    {
        fprintf(stderr, "error: unexpected userdata\n");
        abort();
    }

    flags = JSON_STRICT;
    result = json_unpack_ex(root, &err, flags,
            "{s:o, s:s, s:s}",
            "parameters", &parameters,
            "sequence_a", &sequence_a,
            "sequence_b", &sequence_b);
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
    dp_mat_init(tableau, nrows, ncols);
    sol->mat = tableau;

    char * solution_count_string;
    {
        fmpz_t count;
        fmpz_init(count);
        solve(count, sol, p, A, len_A, B, len_B);
        solution_count_string = fmpz_get_str(NULL, 10, count);
        fmpz_clear(count);
    }

    j_out = json_pack("{s:s}",
            "number_of_optimal_alignments", solution_count_string);

    flint_free(solution_count_string);
    flint_free(A);
    flint_free(B);
    solution_clear(sol);
    model_params_clear(p);
    dp_mat_clear(tableau);

    return j_out;
}



void
solve(fmpz_t res, solution_t sol, const model_params_t p,
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

    tkf91_dp_high(sol, req, mat, expressions_table, generators, A, szA, B, szB);
    count_solutions(res, sol->mat);

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
