/*
 * "Factor Refinement"
 * Eric Bach, James Driscoll, Jeffrey Shallit
 * Computer Sciences Technical Report #883
 * October 1989.
 */
#include "factor_refinement.h"
#include "flint/fmpz_vec.h"


#define fr_node_mref(x) (&(x)->m)

typedef struct fr_node_struct
{
    fmpz m;
    ulong e;
    struct fr_node_struct * next;
} fr_node_struct;
typedef fr_node_struct * fr_node_ptr;


/* functions related to factor refinement nodes */
void fr_node_init(fr_node_ptr x);
void fr_node_init_fmpz_ui(fr_node_ptr x, const fmpz_t m, ulong e);
void fr_node_clear(fr_node_ptr x);
int fr_node_is_one(fr_node_ptr x);
void fr_node_set_fmpz_ui(fr_node_ptr x, const fmpz_t m, ulong e);
void fr_node_get_fmpz_ui(fmpz_t m, ulong * e, fr_node_ptr x);

/* functions related to lists of factor refinement nodes */
slong fr_node_list_length(fr_node_ptr x);
void fr_node_list_pop_front(fr_node_ptr *phead, fr_node_ptr *ptail);
void fr_node_list_concat(fr_node_ptr *phead, fr_node_ptr *ptail,
        fr_node_ptr rhead, fr_node_ptr rtail);
void fr_node_list_clear(fr_node_ptr head);
void fr_node_list_print(fr_node_ptr head);

/* functions related to the actual algorithms of interest */
void pair_refine_unreduced(fr_node_ptr *phead,
        fmpz_t m1, ulong e1, fmpz_t m2, ulong e2);
void remove_ones(fr_node_ptr *phead, fr_node_ptr *ptail, fr_node_ptr ohead);
void pair_refine(fr_node_ptr *phead, fr_node_ptr *ptail,
        fmpz_t m1, ulong e1, fmpz_t m2, ulong e2);
void augment_refinement(fr_node_ptr *phead, fr_node_ptr *ptail,
        const fmpz_t m_jp1, ulong e_jp1,
        fr_node_ptr L_j, fr_node_ptr L_j_tail);

void
fr_node_init(fr_node_ptr x)
{
    fmpz_init(fr_node_mref(x));
    x->e = 0;
    x->next = NULL;
}

void
fr_node_init_fmpz_ui(fr_node_ptr x, const fmpz_t m, ulong e)
{
    fmpz_init_set(fr_node_mref(x), m);
    x->e = e;
    x->next = NULL;
}

void
fr_node_clear(fr_node_ptr x)
{
    fmpz_clear(fr_node_mref(x));
    x->e = 0;
    x->next = NULL;
}

int
fr_node_is_one(fr_node_ptr x)
{
    /* follow the fmpz_pow_ui convention 0^0 = 1 */
    return (x->e == WORD(0) || fmpz_is_one(fr_node_mref(x)));
}

void
fr_node_set_fmpz_ui(fr_node_ptr x, const fmpz_t m, ulong e)
{
    fmpz_set(fr_node_mref(x), m);
    x->e = e;
}

void
fr_node_get_fmpz_ui(fmpz_t m, ulong * e, fr_node_ptr x)
{
    fmpz_set(m, fr_node_mref(x));
    *e = x->e;
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
        *phead = tmp;
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
        *phead = rhead;
        *ptail = rtail;
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
        flint_printf("^%wu ", curr->e);
        curr = curr->next;
    }
    flint_printf("\n");
}

void
pair_refine_unreduced(fr_node_ptr *phead,
        fmpz_t m1, ulong e1,
        fmpz_t m2, ulong e2)
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
    fr_node_init_fmpz_ui(head, m1, e1);

    tail = flint_malloc(sizeof(fr_node_struct));
    fr_node_init_fmpz_ui(tail, m2, e2);

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
                neo->e = curr->e + next->e;

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
        fmpz_t m1, ulong e1,
        fmpz_t m2, ulong e2)
{
    pair_refine_unreduced(phead, m1, e1, m2, e2);
    remove_ones(phead, ptail, *phead);
}

void
augment_refinement(fr_node_ptr *phead, fr_node_ptr *ptail,
        const fmpz_t m_jp1, ulong e_jp1,
        fr_node_ptr L_j, fr_node_ptr L_j_tail)
{
    /* m_jp1 must be positive and greater than 1 */
    /* this function is destructive to the existing refinement L_j */

    fr_node_ptr L_jp1, L_jp1_tail, L_prime, L_prime_tail, neo;
    fmpz_t m;
    ulong e;

    /* initialize (m, e) <- (m_{j+1}, 1), L_{j+1} <- empty list */
    fmpz_init_set(m, m_jp1);
    e = e_jp1;
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
                        fr_node_mref(L_j), L_j->e);

            /* (m, e) <- First(L') */
            fr_node_get_fmpz_ui(m, &e, L_prime);

            /* L' <- Rest(L') */
            fr_node_list_pop_front(&L_prime, &L_prime_tail);

            /* L_{j+1} <- Concat(L_{j+1}, L') */
            fr_node_list_concat(&L_jp1, &L_jp1_tail, L_prime, L_prime_tail);
        }

        /* L_j <- Rest(L_j) */
        fr_node_list_pop_front(&L_j, &L_j_tail);
    }

    /* create a single-element list like [(m, e)] */
    neo = flint_malloc(sizeof(fr_node_struct));
    fr_node_init_fmpz_ui(neo, m, e);

    /* L_{j+1} <- Concat(L_{j+1}, Rest(L_j), [(m, e)]) */
    fr_node_list_pop_front(&L_j, &L_j_tail);
    fr_node_list_concat(&L_jp1, &L_jp1_tail, L_j, L_j_tail);
    fr_node_list_concat(&L_jp1, &L_jp1_tail, neo, neo);

    /* output list of pairs (n_i, e_i) in L_{j+1} with n_i != 1 */
    remove_ones(phead, ptail, L_jp1);

    fmpz_clear(m);
}


int
fmpz_factor_sgn(const fmpz_factor_t f)
{
    int i, s;
    ulong e, neg;
    if (!f->sign)
    {
        return 0;
    }
    neg = f->sign < 0;
    for (i = 0; i < f->num; i++)
    {
        /* follow the fmpz_pow_ui convention 0^0 = 1 */
        e = f->exp[i];
        if (e != WORD(0))
        {
            s = fmpz_sgn(f->p+i);
            if (!s)
            {
                return 0;
            }
            else if (s < 0)
            {
                neg = (neg + e) % 2;
            }
        }
    }
    return neg ? -1 : 1;
}


void
fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f)
{
    int s;
    fr_node_ptr L, L_tail, curr;
    slong i, len;
    ulong e;
    fmpz_t x;

    /* check the sign of f without requiring canonical form */
    s = fmpz_factor_sgn(f);
    if (!s)
    {
        _fmpz_factor_set_length(res, 0);
        res->sign = 0;
        return;
    }

    /* compute the refinement as a linked list */
    fmpz_init(x);
    for (L = L_tail = NULL, i = 0; i < f->num; i++)
    {
        e = f->exp[i];
        if (e != WORD(0))
        {
            fmpz_abs(x, f->p+i);
            if (!fmpz_is_one(x))
            {
                augment_refinement(&L, &L_tail, x, e, L, L_tail);
            }
        }
    }
    fmpz_clear(x);

    /* fill the output with the contents of the linked list */
    len = fr_node_list_length(L);
    _fmpz_factor_fit_length(res, len);
    _fmpz_factor_set_length(res, len);
    res->sign = s;
    for (i = 0, curr = L; i < len; i++, curr = curr->next)
    {
        fmpz_set(res->p+i, fr_node_mref(curr));
        res->exp[i] = curr->e;
    }
    fr_node_list_clear(L);
}
