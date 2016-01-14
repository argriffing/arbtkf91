/*
 * Visualize the tableau.
 */

#include <time.h>

#include "flint/flint.h"
#include "flint/fmpq.h"

#include "jansson.h"

#include "runjson.h"
#include "jsonutil.h"
#include "tkf91_dp_bound.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"
#include "model_params.h"
#include "json_model_params.h"
#include "printutil.h"



void
solve(solution_t sol, int image_mode_full, const char * image_filename,
        const model_params_t p,
        const slong *A, slong len_A, const slong *B, slong len_B);


json_t *run(void * userdata, json_t *root);

json_t *run(void * userdata, json_t *root)
{
    model_params_t p;
    slong len_A, len_B;
    slong *A;
    slong *B;
    solution_t sol;
    int result;
    json_t * parameters;
    const char * sequence_a;
    const char * sequence_b;
    const char * image_filename;
    const char * image_mode;
    int image_mode_full;
    json_error_t err;
    size_t flags;

    if (userdata)
    {
        fprintf(stderr, "error: unexpected userdata\n");
        abort();
    }

    flags = JSON_STRICT;
    json_unpack(root, "{s:o, s:s, s:s, s:s, s:s}",
            "parameters", &parameters,
            "image_mode", &image_mode,
            "image_filename", &image_filename,
            "sequence_a", &sequence_a,
            "sequence_b", &sequence_b);
    if (result)
    {
        fprintf(stderr, "error: on line %d: %s\n", err.line, err.text);
        abort();
    }

    /* read the model parameter values */
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

    /* read the image mode */
    if (strcmp(image_mode, "simple") == 0)
    {
        image_mode_full = 0;
    }
    else if (strcmp(image_mode, "full") == 0)
    {
        image_mode_full = 1;
    }
    else
    {
        fprintf(stderr, "error: expected image_mode : {full,simple}\n");
        abort();
    }

    /* read the two unaligned sequences */

    len_A = strlen(sequence_a);
    A = flint_malloc(len_A * sizeof(slong));
    _fill_sequence_vector(A, sequence_a, len_A);

    len_B = strlen(sequence_b);
    B = flint_malloc(len_B * sizeof(slong));
    _fill_sequence_vector(B, sequence_b, len_B);

    solution_init(sol, len_A + len_B);
    solve(sol, image_mode_full, image_filename, p, A, len_A, B, len_B);

    flint_free(A);
    flint_free(B);
    solution_clear(sol);
    model_params_clear(p);

    return NULL;
}



void
solve(solution_t sol, int image_mode_full, const char * image_filename,
        const model_params_t p,
        const slong *A, slong szA, const slong *B, slong szB)
{
    tkf91_rationals_t r;
    tkf91_expressions_t expressions;
    tkf91_generator_indices_t generators;
    fmpz_mat_t mat;
    expr_ptr * expressions_table;
    request_t req;
    int verbose = 0;
    FILE * file = NULL;
    if (verbose)
    {
        file = stderr;
    }

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
    req->image_mode_full = image_mode_full;
    req->png_filename = image_filename;
    req->trace = 1;

    tkf91_dp_bound(sol, req, mat, expressions_table, generators,
            A, szA, B, szB);

    /* fixme the following block has been moved out of tkf91_dp_bound */
    /* create the tableau png image */
    if (req->png_filename)
    {
        clock_t start = clock();
        if (req->image_mode_full)
        {
            write_tableau_image(
                    req->png_filename, sol->mat, "tkf91 tableau");
        }
        else
        {
            write_simple_tableau_image(
                    req->png_filename, sol->mat, "tkf91 tableau");
        }
        _fprint_elapsed(file, "create tableau png", clock() - start);
    }

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
