/*
 * This executable checks the status of a given alignment
 * in the context of a tkf91 model with given parameter values.
 *
 * Possible outcomes:
 * 1) The alignment is shown to be optimal, both in the sense that no other
 *    alignment has greater probability and that it is the canonical
 *    alignment (according to the 'counter-clockwise traceback' criterion)
 *    among all such alignments.
 * 2) The alignment is shown to have the maximum probability
 *    but is not canonical.
 * 3) The alignment is shown to not have the maximum probability.
 * 4) The status could not be determined.
 *
 * These outcomes are represented by a json object
 * with three members on stdout:
 * {
 * "alignment_is_optimal" : "yes" | "no" | "undetermined",
 * "alignment_is_canonical" : "yes" | "no" | "undetermined",
 * "number_of_optimal_alignments" : <integer-as-string> | "undetermined"
 * }
 * The number of optimal alignments is reported as a json string instead
 * of as a json integer because it is likely to overflow json integer capacity.
 *
 * The json input on stdin is:
 * {
 * "pa_n" : integer, "pa_d" : integer,
 * "pc_n" : integer, "pc_d" : integer,
 * "pg_n" : integer, "pg_d" : integer,
 * "pt_n" : integer, "pt_d" : integer,
 * "lambda_n" : integer, "lambda_d" : integer,
 * "mu_n" : integer, "mu_d" : integer,
 * "tau_n" : integer, "tau_d" : integer,
 * "sequence_a" : string, "sequence_b" : string
 * }
 *
 */

#include "flint/flint.h"
#include "flint/fmpq.h"

#include "jansson.h"

#include "runjson.h"
#include "jsonutil.h"
#include "tkf91_dp_bound.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"
#include "model_params.h"


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


json_t *run(void * userdata, json_t *j_in);

json_t *run(void * userdata, json_t *j_in)
{
    json_t *args;
    json_t *j_out;
    model_params_t p;
    alignment_t aln;
    sequence_pair_t sequences;
    slong len_A, len_B;
    slong *A;
    slong *B;
    solution_t sol;
    breadcrumb_mat_t crumb_mat;

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

    /* model_params_print(p); */

    /* read the two aligned sequences */
    A = _json_object_get_sequence(&len_A, args, "sequence_a");
    B = _json_object_get_sequence(&len_B, args, "sequence_b");

    if (len_A != len_B)
    {
        fprintf(stderr, "error: alignment rows have different lengths\n");
        abort();
    }

    alignment_init(aln, A, B, len_A);
    sequence_pair_init(sequences, aln);

    /* init the tableau mask */
    breadcrumb_mat_init(crumb_mat, sequences->len_A + 1, sequences->len_B + 1);

    /* init the solution object */
    solution_init(sol, sequences->len_A + sequences->len_B);

    /* connect the tableau mask to the solution object */
    sol->pmask = crumb_mat;

    /* do enough of the traceback to get the solution mask */
    /*
    flint_printf("length of sequence A: %wd\n", sequences->len_A);
    flint_printf("length of sequence B: %wd\n", sequences->len_B);
    */
    solve(sol, p, sequences);

    int optimal;
    int canonical;
    breadcrumb_mat_check_alignment(
            &optimal, &canonical, sol->pmask,
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
    breadcrumb_mat_clear(crumb_mat);

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
