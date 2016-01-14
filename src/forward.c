#include "forward.h"
#include "dp.h"


int
dp_forward(dp_mat_t mat, forward_strategy_t strat)
{
    /*
     * Visit only cells where max2 or max3 is interesting.
     * Assume that only the cell of interest and its three neighbors
     * {curr, top, diag, left} are required to update the current cell data.
     * Cell data may be reused.
     * If the visit callback ever returns a nonzero value,
     * then stop the iteration and return that result
     * after clearing the data.
     */
    char *buffer;
    char *curr, *top, *diag, *left;
    char *row, *alt, *tmp;
    size_t sz_cell, active_cell_count;
    int i, j, result;
    slong nrows, ncols;
    dp_t x;

    nrows = dp_mat_nrows(mat);
    ncols = dp_mat_ncols(mat);

    active_cell_count = (size_t) (2 * ncols);
    sz_cell = strat->sz_celldata;

    /* allocate two rows worth of cell data */
    buffer = strat->init(strat->userdata, active_cell_count);
    row = buffer + 0 * ncols * sz_cell;
    alt = buffer + 1 * ncols * sz_cell;

    /* visit interesting cells in row major order */
    result = 0;
    for (i = 0; i < nrows && !result; i++)
    {
        for (j = 0; j < ncols && !result; j++)
        {
            curr = row + j * sz_cell;
            top = diag = left = NULL;
            if (j > 0)
            {
                left = row + (j-1) * sz_cell;
            }
            if (i > 0)
            {
                top = alt + j * sz_cell;
            }
            if (i > 0 && j > 0)
            {
                diag = alt + (j-1) * sz_cell;
            }
            x = *dp_mat_entry(mat, i, j);
            if (x & (DP_MAX2 | DP_MAX3))
            {
                result = strat->visit(strat->userdata,
                        mat, i, j, curr, top, diag, left);
            }
        }
        tmp = row;
        row = alt;
        alt = tmp;
    }

    strat->clear(strat->userdata, buffer, active_cell_count);
    return result;
}
