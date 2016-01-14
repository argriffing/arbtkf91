/*
 * This executable checks the status of a given alignment
 * in the context of a tkf91 model with given parameter values.
 *
 * Output format:
 * {
 * "alignment_is_optimal" : "yes" | "no" | "undetermined",
 * "alignment_is_canonical" : "yes" | "no" | "undetermined",
 * "number_of_optimal_alignments" : <integer-as-string> | "undetermined"
 * }
 * The number of optimal alignments is reported as a json string instead
 * of as a json integer because it is likely to overflow json integer capacity.
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
#include "bound_mat.h"
#include "count_solutions.h"
#include "printutil.h"


typedef struct
{
    slong *A;
    slong *B;
    slong len;
} alignment_struct;
typedef alignment_struct alignment_t[1];

static void alignment_init(alignment_t x, slong *A, slong *B, slong len);
static void alignment_clear(alignment_t x);

void
alignment_init(alignment_t x, slong *A, slong *B, slong len)
{
    x->A = A;
    x->B = B;
    x->len = len;
}

void
alignment_clear(alignment_t x)
{
    flint_free(x->A);
    flint_free(x->B);
}



typedef struct
{
    slong *A;
    slong *B;
    slong len_A;
    slong len_B;
} sequence_pair_struct;
typedef sequence_pair_struct sequence_pair_t[1];

static void sequence_pair_init(sequence_pair_t x, alignment_t aln);
static void sequence_pair_clear(sequence_pair_t x);

void
sequence_pair_init(sequence_pair_t x, alignment_t aln)
{
    slong i;
    slong value;
    x->A = malloc(aln->len * sizeof(slong));
    x->B = malloc(aln->len * sizeof(slong));
    x->len_A = 0;
    x->len_B = 0;
    
    for (i = 0; i < aln->len; i++)
    {
        if (aln->A[i] >= 0)
        {
            value = aln->A[i];
            /* flint_printf("A[%wd] : %wd\n", x->len_A, value); */
            x->A[(x->len_A)++] = value;
        }
        if (aln->B[i] >= 0)
        {
            value = aln->B[i];
            /* flint_printf("B[%wd] : %wd\n", x->len_B, value); */
            x->B[(x->len_B)++] = aln->B[i];
        }
    }
}

void
sequence_pair_clear(sequence_pair_t x)
{
    free(x->A);
    free(x->B);
}





void
solve(solution_t sol, const model_params_t p,
        const sequence_pair_t sequences);


json_t *run(void * userdata, json_t *root);

json_t *run(void * userdata, json_t *root)
{
    json_t *j_out;
    model_params_t p;
    alignment_t aln;
    sequence_pair_t sequences;
    slong len_A, len_B;
    slong *A;
    slong *B;
    solution_t sol;
    dp_mat_t tableau;
    int result;
    const char * sequence_a;
    const char * sequence_b;
    json_t * parameters;
    json_error_t err;
    size_t flags;
    slong nrows, ncols;

    if (userdata)
    {
        fprintf(stderr, "error: unexpected userdata\n");
        abort();
    }

    flags = JSON_STRICT;
    result = json_unpack_ex(root, &err, flags, "{s:o, s:s, s:s}",
            "parameters", &parameters,
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

    /* read the two aligned sequences */

    len_A = strlen(sequence_a);
    A = flint_malloc(len_A * sizeof(slong));
    _fill_sequence_vector(A, sequence_a, len_A);

    len_B = strlen(sequence_b);
    B = flint_malloc(len_B * sizeof(slong));
    _fill_sequence_vector(B, sequence_b, len_B);

    if (len_A != len_B)
    {
        fprintf(stderr, "error: alignment rows have different lengths\n");
        abort();
    }

    alignment_init(aln, A, B, len_A);
    sequence_pair_init(sequences, aln);

    /* init the solution object */
    nrows = sequences->len_A + 1;
    ncols = sequences->len_B + 1;
    solution_init(sol, nrows + ncols);

    /* init the tableau and connect it to the solution object */
    dp_mat_init(tableau, nrows, ncols);
    sol->mat = tableau;

    /* do enough of the traceback to get the solution mask */
    solve(sol, p, sequences);

    int optimal;
    int canonical;
    dp_mat_check_alignment(
            &optimal, &canonical, sol->mat,
            aln->A, aln->B, aln->len);

    char * solution_count_string;
    char * _solution_count_string = NULL;
    char undetermined[] = "undetermined";
    solution_count_string = undetermined;
    if (sol->has_best_tie_count)
    {
        _solution_count_string = fmpz_get_str(NULL, 10, sol->best_tie_count);
        solution_count_string = _solution_count_string;
    }

    j_out = json_pack("{s:s, s:s, s:s}",
            "alignment_is_optimal", (optimal ? "yes" : "no"),
            "alignment_is_canonical", (canonical ? "yes" : "no"),
            "number_of_optimal_alignments", solution_count_string);

    solution_clear(sol);
    model_params_clear(p);
    alignment_clear(aln);
    sequence_pair_clear(sequences);
    flint_free(_solution_count_string);
    dp_mat_clear(tableau);

    return j_out;
}



void
solve(solution_t sol, const model_params_t p,
        const sequence_pair_t sequences)
{
    tkf91_rationals_t r;
    tkf91_expressions_t expressions;
    tkf91_generator_indices_t generators;
    fmpz_mat_t mat;
    expr_ptr * expressions_table;
    request_t req;
    slong *A = sequences->A;
    slong szA = sequences->len_A;
    slong *B = sequences->B;
    slong szB = sequences->len_B;
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
    req->png_filename = NULL;
    req->trace = 1;

    tkf91_dp_bound(
            sol, req, mat, expressions_table, generators,
            A, szA, B, szB);

    /* fixme the following code block has been moved from tkf91_dp_bound */
    /* symbolic verification */
    /* this should be updated to use the forward pass framework */
    {
        clock_t start = clock();
        tkf91_dp_verify_symbolically(
                &sol->optimality_flag, 
                mat, generators, sol->mat,
                expressions_table,
                A, B);
        _fprint_elapsed(file, "symbolic verification", clock() - start);
    }

    /* fixme the following code block has been moved from tkf91_dp_bound */
    /* count the solutions */
    {
        clock_t start = clock();
        fmpz_t solution_count;
        fmpz_init(solution_count);

        count_solutions(solution_count, sol->mat);

        fmpz_set(sol->best_tie_count, solution_count);
        sol->has_best_tie_count = 1;

        fmpz_clear(solution_count);
        _fprint_elapsed(file, "counting solutions", clock() - start);
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
