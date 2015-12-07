/*
 * "Factor Refinement"
 * Eric Bach, James Driscoll, Jeffrey Shallit
 * Computer Sciences Technical Report #883
 * October 1989.
 *
 * GPL implementation by Alex G. 2015.
 */

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"


/* factor refinement node */

#define fr_node_mref(x) (&(x)->m)
#define fr_node_eref(x) (&(x)->e)

typedef struct fr_node_struct
{
    fmpz m;
    fmpz e;
    struct fr_node_struct * next;
} fr_node_struct;

typedef fr_node_struct * fr_node_ptr;

void
fr_node_init(fr_node_ptr x)
{
    fmpz_init(fr_node_mref(x));
    fmpz_init(fr_node_eref(x));
    x->next = NULL;
}

void
fr_node_init_fmpz_fmpz(fr_node_ptr x, const fmpz_t m, const fmpz_t e)
{
    fmpz_init_set(fr_node_mref(x), m);
    fmpz_init_set(fr_node_eref(x), e);
    x->next = NULL;
}

void
fr_node_clear(fr_node_ptr x)
{
    fmpz_clear(fr_node_mref(x));
    fmpz_clear(fr_node_eref(x));
    x->next = NULL;
}

int
fr_node_is_one(fr_node_ptr x)
{
    return fmpz_is_one(fr_node_mref(x));
}

void
fr_node_set_fmpz_fmpz(fr_node_ptr x, const fmpz_t m, const fmpz_t e)
{
    fmpz_set(fr_node_mref(x), m);
    fmpz_set(fr_node_eref(x), e);
}

void
fr_node_get_fmpz_fmpz(fmpz_t m, fmpz_t e, fr_node_ptr x)
{
    fmpz_set(m, fr_node_mref(x));
    fmpz_set(e, fr_node_eref(x));
}

slong
fr_node_list_length(fr_node_ptr x)
{
    slong count;
    count = 0;
    while (x)
    {
        count++;
        x = x->next;
    }
    return count;
}

void
fr_node_list_pop_front(fr_node_ptr *phead, fr_node_ptr *ptail)
{
    fr_node_ptr tmp;

    if (phead == ptail)
    {
        flint_printf("aliasing issue...\n");
        abort();
    }

    if (*phead)
    {

        if (*phead == *ptail)
        {
            *ptail = NULL;
        }

        tmp = (*phead)->next;
        fr_node_clear(*phead);
        flint_free(*phead);
        (*phead) = tmp;
    }
}

void
fr_node_list_concat(fr_node_ptr *phead, fr_node_ptr *ptail,
        fr_node_ptr rhead, fr_node_ptr rtail)
{
    /* if the second list is empty then do nothing */
    if (!rhead)
    {
        return;
    }

    /* if the first list is empty then use the second list */
    if (!(*phead))
    {
        (*phead) = rhead;
        (*ptail) = rtail;
        return;
    }

    /* neither list is empty so connect the tail to the head */
    (*ptail)->next = rhead;
    *ptail = rtail;
}

void
fr_node_list_clear(fr_node_ptr head)
{
    fr_node_ptr curr, next;
    curr = head;
    while (curr)
    {
        next = curr->next;
        fr_node_clear(curr);
        flint_free(curr);
        curr = next;
    }
}

void
fr_node_list_print(fr_node_ptr head)
{
    fr_node_ptr curr;

    curr = head;
    while (curr)
    {
        fmpz_print(fr_node_mref(curr));
        flint_printf("^");
        fmpz_print(fr_node_eref(curr));
        flint_printf(" ");
        curr = curr->next;
    }
    flint_printf("\n");
}

void
pair_refine_unreduced(fr_node_ptr *phead,
        fmpz_t m1, fmpz_t e1,
        fmpz_t m2, fmpz_t e2)
{
    fr_node_ptr head, tail, curr, next, neo;
    fmpz_t d;
    int boring;

    if (fmpz_is_one(m1) && fmpz_is_one(m2))
    {
        *phead = NULL;
        return;
    }

    fmpz_init(d);

    head = flint_malloc(sizeof(fr_node_struct));
    fr_node_init_fmpz_fmpz(head, m1, e1);

    tail = flint_malloc(sizeof(fr_node_struct));
    fr_node_init_fmpz_fmpz(tail, m2, e2);

    head->next = tail;

    boring = 0;
    while (!boring)
    {
        curr = head;
        next = curr->next;

        boring = 1;
        while (next)
        {
            if (!fr_node_is_one(curr) && !fr_node_is_one(next))
            {
                fmpz_gcd(d, fr_node_mref(curr), fr_node_mref(next));
                fmpz_divexact(fr_node_mref(curr), fr_node_mref(curr), d);
                fmpz_divexact(fr_node_mref(next), fr_node_mref(next), d);

                neo = flint_malloc(sizeof(fr_node_struct));
                fr_node_init(neo);
                fmpz_set(fr_node_mref(neo), d);
                fmpz_add(fr_node_eref(neo),
                         fr_node_eref(curr), fr_node_eref(next));

                curr->next = neo;
                neo->next = next;

                next = neo;
                boring = 0;
            }
            else
            {
                curr = next;
                next = next->next;
            }
        }
    }

    fmpz_clear(d);
    *phead = head;
}

void
remove_ones(fr_node_ptr *phead, fr_node_ptr *ptail, fr_node_ptr ohead)
{
    fr_node_ptr head, curr, ocurr, onext;

    if (!ohead)
    {
        *phead = NULL;
        *ptail = NULL;
        return;
    }

    head = NULL;
    curr = NULL;
    ocurr = ohead;
    while (ocurr)
    {
        onext = ocurr->next;

        if (fr_node_is_one(ocurr))
        {
            fr_node_clear(ocurr);
            flint_free(ocurr);
        }
        else
        {
            if (!head)
            {
                head = ocurr;
                curr = ocurr;
            }
            else
            {
                curr->next = ocurr;
                curr = curr->next;
            }
        }
        ocurr = onext;
    }
    curr->next = NULL;

    *phead = head;
    *ptail = curr;
}

void
pair_refine(fr_node_ptr *phead, fr_node_ptr *ptail,
        fmpz_t m1, fmpz_t e1,
        fmpz_t m2, fmpz_t e2)
{
    pair_refine_unreduced(phead, m1, e1, m2, e2);
    remove_ones(phead, ptail, *phead);
}

void
pair_refine_si(fr_node_ptr *phead,
        slong m1, slong e1,
        slong m2, slong e2)
{
    fmpz_t zm1, ze1, zm2, ze2;
    fr_node_ptr tail;

    fmpz_init_set_si(zm1, m1);
    fmpz_init_set_si(ze1, e1);
    fmpz_init_set_si(zm2, m2);
    fmpz_init_set_si(ze2, e2);

    pair_refine(phead, &tail, zm1, ze1, zm2, ze2);

    fmpz_clear(zm1);
    fmpz_clear(ze1);
    fmpz_clear(zm2);
    fmpz_clear(ze2);
}


void
augment_refinement(fr_node_ptr *phead, fr_node_ptr *ptail,
        const fmpz_t m_jp1,
        fr_node_ptr L_j, fr_node_ptr L_j_tail)
{
    /* this function is destructive to the existing refinement L_j */

    fr_node_ptr L_jp1, L_jp1_tail, L_prime, L_prime_tail, neo;
    fmpz_t m, e;

    /* initialize (m, e) <- (m_{j+1}, 1), L_{j+1} <- empty list */
    fmpz_init_set(m, m_jp1);
    fmpz_init_set_ui(e, 1);
    L_jp1 = NULL;
    L_jp1_tail = NULL;
    L_prime = NULL;
    L_prime_tail = NULL;

    while (L_j && !fmpz_is_one(m))
    {
        if (!fr_node_is_one(L_j))
        {
            /* L' <- Pair-Refine((m, e), First(L_j)) */
            pair_refine(&L_prime, &L_prime_tail, m, e,
                        fr_node_mref(L_j), fr_node_eref(L_j));

            /* (m, e) <- First(L') */
            fr_node_get_fmpz_fmpz(m, e, L_prime);

            /* L' <- Rest(L') */
            fr_node_list_pop_front(&L_prime, &L_prime_tail);

            /* L_{j+1} <- Concat(L_{j+1}, L') */
            fr_node_list_concat(&L_jp1, &L_jp1_tail, L_prime, L_prime_tail);
        }

        /* L_j <- Rest(L_j) */
        fr_node_list_pop_front(&L_j, &L_j_tail);
    }

    /* create a single-element list like (m, e) */
    neo = flint_malloc(sizeof(fr_node_struct));
    fr_node_init_fmpz_fmpz(neo, m, e);

    /* L_{j+1} <- Concat(L_{j+1}, Rest(L_j), (m, e)) */
    fr_node_list_pop_front(&L_j, &L_j_tail);
    fr_node_list_concat(&L_jp1, &L_jp1_tail, L_j, L_j_tail);
    fr_node_list_concat(&L_jp1, &L_jp1_tail, neo, neo);

    /* output list of pairs (n_i, e_i) in L_{j+1} with n_i != 1 */
    remove_ones(phead, ptail, L_jp1);

    fmpz_clear(m);
    fmpz_clear(e);
}

void
factor_refine(fmpz **ybase, fmpz **yexp, slong *ylen,
              const fmpz *x, const slong xlen)
{
    fr_node_ptr L, L_tail, curr;
    slong i;

    /* compute the refinement as a linked list */
    L = NULL;
    L_tail = NULL;
    for (i = 0; i < xlen; i++)
    {
        if (!fmpz_is_one(x+i))
        {
            augment_refinement(&L, &L_tail, x+i, L, L_tail);
        }
    }

    /* convert the linked list to a more standard format */
    *ylen = fr_node_list_length(L);
    *ybase = _fmpz_vec_init(*ylen);
    *yexp = _fmpz_vec_init(*ylen);
    curr = L;
    i = 0;
    while (curr)
    {
        fmpz_set((*ybase)+i, fr_node_mref(curr));
        fmpz_set((*yexp)+i, fr_node_eref(curr));
        curr = curr->next;
        i++;
    }

    fr_node_list_clear(L);
}

void
_fmpz_vec_randtest_pos(fmpz * f, flint_rand_t state,
        slong len, mp_bitcnt_t bits)
{
    slong i;
    _fmpz_vec_randtest_unsigned(f, state, len, bits-1);
    for (i = 0; i < len; i++)
    {
        fmpz_add_ui(f+i, f+i, 1);
    }
}

int test_factor_refinement()
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("factor_refinement....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        fmpz *x, *ybase, *yexp;
        slong xlen, ylen;

        mp_bitcnt_t bits;

        x = ybase = yexp = NULL;
        xlen = ylen = 0;

        bits = n_randint(state, 100) + 2;
        xlen = n_randint(state, 100) + 1;

        /* create the random input vector of positive integers */
        x = _fmpz_vec_init(xlen);
        _fmpz_vec_randtest_pos(x, state, xlen, bits);

        factor_refine(&ybase, &yexp, &ylen, x, xlen);

        /* check that products are equal */
        {
            fmpz_t a, b, p;
            slong j, u;

            fmpz_init_set_ui(a, 1);
            for (j = 0; j < xlen; j++)
            {
                fmpz_mul(a, a, x+j);
            }

            fmpz_init_set_ui(b, 1);
            fmpz_init(p);
            for (j = 0; j < ylen; j++)
            {
                u = fmpz_get_ui(yexp+j);
                fmpz_pow_ui(p, ybase+j, u);
                fmpz_mul(b, b, p);
            }

            if (!fmpz_equal(a, b))
            {
                flint_printf("FAIL:\n");
                fmpz_print(a); flint_printf(" ");
                fmpz_print(b); flint_printf("\n");
                abort();
            }

            fmpz_clear(a);
            fmpz_clear(b);
            fmpz_clear(p);
        }

        /* check that elements of the base are pairwise coprime */
        {
            slong u, v;
            fmpz_t g;
            fmpz_init(g);
            for (u = 0; u < ylen; u++)
            {
                for (v = 0; v < u; v++)
                {
                    fmpz_gcd(g, ybase+u, ybase+v);
                    if (!fmpz_is_one(g))
                    {
                        flint_printf("FAIL:\n");
                        fmpz_print(ybase+u); flint_printf(" ");
                        fmpz_print(ybase+v); flint_printf("\n");
                        abort();
                    }
                }
            }
            fmpz_clear(g);
        }

        /* check that each input is a product of powers of outputs */
        {
            slong u, v;
            fmpz_t a;
            fmpz_init(a);
            for (u = 0; u < xlen; u++)
            {
                fmpz_set(a, x+u);
                for (v = 0; v < ylen && !fmpz_is_one(a); v++)
                {
                    while (fmpz_divisible(a, ybase+v))
                    {
                        fmpz_divexact(a, a, ybase+v);
                    }
                }
                if (!fmpz_is_one(a))
                {
                    flint_printf("FAIL:\n");
                    fmpz_print(ybase+u); flint_printf("\n");
                    abort();
                }
            }
            fmpz_clear(a);
        }

        _fmpz_vec_clear(x, xlen);
        _fmpz_vec_clear(ybase, ylen);
        _fmpz_vec_clear(yexp, ylen);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}


int main(int argc, char *argv[])
{
    /*
    fr_node_ptr head;
    pair_refine_si(&head, 3*4*5, 2, 5*3*8*11, 3);

    fr_node_list_print(head);
    fr_node_list_clear(head);
    */

    /*
    int success;
    
    fmpz *x, *ybase, *yexp;
    slong xlen, ylen;

    x = ybase = yexp = NULL;
    xlen = ylen = 0;

    flint_printf("enter a vector:\n");
    success = _fmpz_vec_read(&x, &xlen);
    if (!success)
    {
        flint_printf("failed to read the vector\n");
    }

    flint_printf("you entered this vector:\n");
    _fmpz_vec_print(x, xlen);
    flint_printf("\n");

    flint_printf("refining the factorization...\n");
    factor_refine(&ybase, &yexp, &ylen, x, xlen);

    flint_printf("base: ");
    _fmpz_vec_print(ybase, ylen);
    flint_printf("exp: ");
    _fmpz_vec_print(yexp, ylen);

    _fmpz_vec_clear(x, xlen);
    _fmpz_vec_clear(ybase, ylen);
    _fmpz_vec_clear(yexp, ylen);

    return 0;
    */
    
    return test_factor_refinement();
}
