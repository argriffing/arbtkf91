#ifndef TKF91_GENERATOR_VECS_H
#define TKF91_GENERATOR_VECS_H

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#include "tkf91_generator_indices.h"


/*
 * Vec-ification of tkf91 generators.
 * The members of this structure will point to rows of an fmpz matrix.
 * If my ideas are right, the fmpz matrix will be the inverse of the
 * transform matrix associated with the hermite normal form transformation
 * of the generator matrix.
 * Each of the vectors will have length equal to the rank of the Hermite 
 * normal form, but the 'hermitification' procedures do not need to
 * know or care about this length.
 */
typedef struct
{
    fmpz * m1_00;
    fmpz * m0_10;
    fmpz * m0_i0_incr[4];
    fmpz * m2_01;
    fmpz * m2_0j_incr[4];
    fmpz * c0_incr[4];
    fmpz * c1_incr[16];
    fmpz * c2_incr[4];
    fmpz_mat_t M;
} tkf91_generator_vecs_struct;
typedef tkf91_generator_vecs_struct tkf91_generator_vecs_t[1];



#ifdef __cplusplus
extern "C" {
#endif


void _fmpz_mat_hnf_inverse_transform(
        fmpz_mat_t H, fmpz_mat_t V, slong * prank, const fmpz_mat_t A);


static __inline__ slong
tkf91_generator_vecs_rank(tkf91_generator_vecs_t h)
{
    return fmpz_mat_ncols(h->M);
}


void tkf91_generator_vecs_init(
        tkf91_generator_vecs_t h,
        const tkf91_generator_indices_t g,
        const fmpz_mat_t V,
        slong rank);

void tkf91_generator_vecs_clear(tkf91_generator_vecs_t h);


#ifdef __cplusplus
}
#endif

#endif
