/*
 * Given a list of rational numbers, define some kind of decomposition.
 * Compute a list of co-prime integers.
 * Write each of the provided rational numbers as products of integer powers
 * of the co-prime integers in the computed list.
 *
 * Reference:
 * Factor Refinement
 * Eric Bach, James Driscoll, Jeffrey Shallit
 * Computer Sciences Technical Report #883
 * October 1989.
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
    /* if the first list is empty then use the second list */
    if (!(*phead))
    {
        (*phead) = rhead;
        (*ptail) = rtail;
        return;
    }
    
    /* if the second list is empty then use the first list */
    if (!rhead)
    {
        return;
    }

    /* if neither list is empty then connect the tail to the head */
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
    /* phead will be set to the head of a list without ones */

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
                curr = head;
            }
            else
            {
                curr->next = ocurr;
                curr = ocurr;
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

    flint_printf("before removing ones in pair refinement...\n");
    fr_node_list_print(*phead);
    flint_printf("\n");

    remove_ones(phead, ptail, *phead);

    flint_printf("after removing ones in pair refinement...\n");
    fr_node_list_print(*phead);
    flint_printf("\n");
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

            /* FIXME this is for debugging ... */
            fr_node_list_print(L_prime);

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
    fr_node_list_pop_front(&L_j, &L_j_tail); /* FIXME not sure about this... */

    flint_printf("before first concatenation: ");
    fr_node_list_print(L_jp1);

    fr_node_list_concat(&L_jp1, &L_jp1_tail, L_j, L_j_tail);

    flint_printf("before second concatenation: ");
    fr_node_list_print(L_jp1);

    fr_node_list_concat(&L_jp1, &L_jp1_tail, neo, neo);

    flint_printf("after second concatenation: ");
    fr_node_list_print(L_jp1);

    /* output list of pairs (n_i, e_i) in L_{j+1} with n_i != 1 */
    remove_ones(phead, ptail, L_jp1);

    flint_printf("after removing ones in augment_refinement: ");
    fr_node_list_print(*phead);

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
        augment_refinement(&L, &L_tail, x+i, L, L_tail);
    }

    flint_printf("linked list before converting to vec: \n");
    fr_node_list_print(L);

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

int main(int argc, char *argv[])
{

    /*
    fr_node_ptr head;
    pair_refine_si(&head, 3*4*5, 2, 5*3*8*11, 3);

    fr_node_list_print(head);
    fr_node_list_clear(head);
    */

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
}
