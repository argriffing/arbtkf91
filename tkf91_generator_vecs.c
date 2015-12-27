#include "flint/fmpz"
#include "flint/fmpz_mat"

#include "tkf91_generator_indices.h"
#include "tkf91_generator_vecs.h"

fmpz * _vecify(slong i, const fmpz_mat_t M);

fmpz *
_vecify(slong i, const fmpz_mat_t M)
{
    return M->rows[i];
}

void
vecify_tkf91_generators(
        tkf91_hermite_generators_t h,
        const tkf91_generator_indices_t g,
        const fmpz_mat_t M)
{
    slong i, j;
    h->m1_00 = _vecify(g->m1_00, M);
    h->m0_10 = _vecify(g->m0_10, M);
    h->m2_01 = _vecify(g->m2_01, M);
    for (i = 0; i < 4; i++)
    {
        h->m0_i0_incr[i] = _vecify(g->m0_i0_incr[i], M);
        h->m2_0j_incr[i] = _vecify(g->m2_0j_incr[i], M);
        h->c0_incr[i] = _vecify(g->c0_incr[i], M);
        for (j = 0; j < 4; j++)
        {
            h->c1_incr[i*4+j] = _vecify(g->c1_incr[i*4+j], M);
        }
        h->c2_incr[i] = _vecify(g->c2_incr[i], M);
    }
}
