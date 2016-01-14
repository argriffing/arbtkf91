#ifndef TKF91_GENERATOR_VECS_H
#define TKF91_GENERATOR_VECS_H

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#include "arb.h"

#include "expressions.h"
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
typedef tkf91_generator_vecs_struct * tkf91_generator_vecs_ptr;



#ifdef __cplusplus
extern "C" {
#endif


void
_arb_vec_dot_fmpz_vec(
        arb_t res, arb_srcptr vec1, const fmpz * vec2, slong len2, slong prec);


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

/*
 * Compute log probabilities corresponding to linear integer combinations
 * of the basis expressions.
 * This function could be used for verification stages of alignment algorithms
 * that are based on (mag_t low, mag_t high) or arb_t intervals.
 *
 * H : the Hermite normal form of G, where G is
 *     the (#generators x #expressions) matrix of integer exponents.
 * expressions_table : array of basis expression objects
 * rank : rank of G and H
 * level : log2 of the bits of precision for intermediate calculations
 */
void compute_hlogy(arb_ptr res, const fmpz_mat_t H,
        expr_ptr * expressions_table, slong rank, slong level);


#ifdef __cplusplus
}
#endif

#endif
