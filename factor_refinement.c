/*
 * "Factor Refinement"
 * Eric Bach, James Driscoll, Jeffrey Shallit
 * Computer Sciences Technical Report #883
 * October 1989.
 *
 * GPL implementation by Alex G. 2015.
 */

#include "factor_refinement.h"
#include "flint/fmpz_vec.h"


#define fr_node_mref(x) (&(x)->m)
#define fr_node_eref(x) (&(x)->e)

typedef struct fr_node_struct
{
    fmpz m;
    fmpz e;
    struct fr_node_struct * next;
} fr_node_struct;

typedef fr_node_struct * fr_node_ptr;


/* the harsher gcc warning flags are making me add some declarations */

/* functions directly related to factor refinement nodes */
void fr_node_init(fr_node_ptr x);
void fr_node_init_fmpz_fmpz(fr_node_ptr x, const fmpz_t m, const fmpz_t e);
void fr_node_clear(fr_node_ptr x);
int fr_node_is_one(fr_node_ptr x);
void fr_node_set_fmpz_fmpz(fr_node_ptr x, const fmpz_t m, const fmpz_t e);
void fr_node_get_fmpz_fmpz(fmpz_t m, fmpz_t e, fr_node_ptr x);

/* functions related to lists of factor refinement nodes */
slong fr_node_list_length(fr_node_ptr x);
void fr_node_list_pop_front(fr_node_ptr *phead, fr_node_ptr *ptail);
void fr_node_list_concat(fr_node_ptr *phead, fr_node_ptr *ptail,
        fr_node_ptr rhead, fr_node_ptr rtail);
void fr_node_list_clear(fr_node_ptr head);
void fr_node_list_print(fr_node_ptr head);

/* functions related to the actual algorithms of interest */
void pair_refine_unreduced(fr_node_ptr *phead,
        fmpz_t m1, fmpz_t e1, fmpz_t m2, fmpz_t e2);
void remove_ones(fr_node_ptr *phead, fr_node_ptr *ptail, fr_node_ptr ohead);
void pair_refine(fr_node_ptr *phead, fr_node_ptr *ptail,
        fmpz_t m1, fmpz_t e1, fmpz_t m2, fmpz_t e2);
void augment_refinement(fr_node_ptr *phead, fr_node_ptr *ptail,
        const fmpz_t m_jp1,
        fr_node_ptr L_j, fr_node_ptr L_j_tail);

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
factor_refinement(fmpz **ybase, fmpz **yexp, slong *ylen,
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
