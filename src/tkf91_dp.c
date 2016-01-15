#include <stdlib.h>

#include "flint/flint.h"
#include "flint/fmpz.h"

#include "arb.h"

#include "tkf91_dp.h"

void
solution_init(solution_t x, slong aln_maxlen)
{
    x->A = flint_calloc(aln_maxlen + 1, sizeof(char));
    x->B = flint_calloc(aln_maxlen + 1, sizeof(char));
    x->len = -1;
    arb_init(x->log_probability);
    x->optimality_flag = 0;
    x->mat = NULL;
}

void
solution_clear(solution_t x)
{
    flint_free(x->A);
    flint_free(x->B);
    arb_clear(x->log_probability);
}

void
solution_fprint(FILE * file, const solution_t x)
{
    /* print the two sequences */
    flint_fprintf(file, "alignment:\n");
    if (x->A)
    {
        flint_fprintf(file, "%s\n", x->A);
    }
    if (x->B)
    {
        flint_fprintf(file, "%s\n", x->B);
    }

    /* report the log probability of the best alignment */
    flint_fprintf(file, "alignment log probability : ");
    arb_fprintd(file, x->log_probability, 15);
    flint_fprintf(file, "\n");
    flint_fprintf(file, "... and in a more nerdy format : ");
    arb_fprint(file, x->log_probability);
    flint_fprintf(file, "\n");
}
