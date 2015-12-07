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
fr_node_clear(fr_node_ptr x)
{
    fmpz_clear(fr_node_mref(x));
    fmpz_clear(fr_node_eref(x));
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
fr_node_clear_list(fr_node_ptr head)
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
pair_refine(fr_node_ptr *phead,
        fmpz_t m1, fmpz_t e1,
        fmpz_t m2, fmpz_t e2)
{
    fr_node_ptr head, tail, curr, next, neo;
    fmpz_t d;
    int boring;

    if (fmpz_is_one(m1) && fmpz_is_one(m2))
    {
        *phead = NULL;
    }

    fmpz_init(d);

    head = flint_malloc(sizeof(fr_node_struct));
    fr_node_init(head);
    fr_node_set_fmpz_fmpz(head, m1, e1);

    tail = flint_malloc(sizeof(fr_node_struct));
    fr_node_init(tail);
    fr_node_set_fmpz_fmpz(tail, m2, e2);

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
pair_refine_si(fr_node_ptr *phead,
        slong m1, slong e1,
        slong m2, slong e2)
{
    fmpz_t zm1, ze1, zm2, ze2;

    fmpz_init_set_si(zm1, m1);
    fmpz_init_set_si(ze1, e1);
    fmpz_init_set_si(zm2, m2);
    fmpz_init_set_si(ze2, e2);

    pair_refine(phead, zm1, ze1, zm2, ze2);

    fmpz_clear(zm1);
    fmpz_clear(ze1);
    fmpz_clear(zm2);
    fmpz_clear(ze2);
}

fr_node_print_list(fr_node_ptr head)
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

int main(int argc, char *argv[])
{
    fr_node_ptr head;
    pair_refine_si(&head, 3*4*5, 2, 5*3*8*11, 3);

    fr_node_print_list(head);

    return 0;
}

