/*
 * Visualize a dynamic programming tableau
 * by creating a .png image with libpng.
 *
 * sudo apt-get install libpng-dev
 * link to -lpng
 * the header is png.h
 */

#include <stdio.h>

#include "flint/flint.h"
#include "png.h"

#include "breadcrumbs.h"
#include "vis.h"

typedef struct
{
    png_bytep data;
    slong r;
    slong c;
    slong sz_data;
} buf_struct;
typedef buf_struct buf_t[1];

static void buf_init(buf_t buf, slong nrows, slong ncols);
static void buf_clear(buf_t buf);
static void buf_zero(buf_t buf);
static png_bytep buf_pixel_row(buf_t buf, slong di);
static png_bytep buf_pixel(buf_t buf, slong j, slong di, slong dj);
static void buf_set_rgb(buf_t buf, slong j, slong di, slong dj,
        png_byte r, png_byte g, png_byte b);

void
buf_init(buf_t buf, slong nrows, slong ncols)
{
    buf->sz_data = 5 * ncols * 5 * 3 * sizeof(png_byte);
    buf->data = malloc(buf->sz_data);
    buf->r = nrows;
    buf->c = ncols;
}

void
buf_clear(buf_t buf)
{
    free(buf->data);
}

png_bytep
buf_pixel_row(buf_t buf, slong di)
{
    return buf->data + di * (buf->c * 5 * 3);
}

png_bytep
buf_pixel(buf_t buf, slong j, slong di, slong dj)
{
    /* the buffer has rgb values for each 5 pixel x 5 pixel cell
     * in one row of the tableau.
     * j is the column of the dp cell
     * di and dj index the row and column of the pixel within the dp cell */
    return buf_pixel_row(buf, di) + j * (5 * 3) + dj * 3;
}

void
buf_set_rgb(buf_t buf, slong j, slong di, slong dj,
        png_byte r, png_byte g, png_byte b)
{
    png_bytep pixel;
    pixel = buf_pixel(buf, j, di, dj);
    pixel[0] = r;
    pixel[1] = g;
    pixel[2] = b;
}

void
buf_zero(buf_t buf)
{
    memset(buf->data, 0, buf->sz_data);
}



int write_tableau_image(
        char * filename, const breadcrumb_mat_t mat, char * title)
{
    FILE *fout;
    png_structp png_ptr;
    png_infop info_ptr;
    int code, width, height;
    png_byte r, g, b;
    breadcrumb_t curr, top, diag, left;

    fout = NULL;
    png_ptr = NULL;
    info_ptr = NULL;
    code = 0;

    slong nrows, ncols;

    nrows = breadcrumb_mat_nrows(mat);
    ncols = breadcrumb_mat_ncols(mat);

    width = (int) ncols * 5 - 2;
    height = (int) nrows * 5 - 2;

    fout = fopen(filename, "wb");
    if (fout == NULL)
    {
        fprintf(stderr, "failed to open %s for writing\n", filename);
        code = 1;
        goto end;
    }
    
    png_ptr = png_create_write_struct(
            PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL)
    {
        fprintf(stderr, "could not allocate write struct\n");
        code = 1;
        goto end;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL)
    {
        fprintf(stderr, "could not allocate info struct\n");
        code = 1;
        goto end;
    }

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        fprintf(stderr, "error during png creation\n");
        code = 1;
        goto end;
    }

    png_init_io(png_ptr, fout);

    /* write header (8 bit color depth) */
    png_set_IHDR(png_ptr, info_ptr, width, height,
            8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != NULL)
    {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = "Title";
        title_text.text = title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    /* write image data one row at a time */
    breadcrumb_t want_any = (CRUMB_WANT2 | CRUMB_WANT3 | CRUMB_CONTENDER);
    int i, j, di, dj;
    buf_t buf;
    buf_init(buf, nrows, ncols);
    for (i = 0; i < nrows; i++)
    {
        /* reset all entries of the pixel buffer to zero */
        buf_zero(buf);

        for (j = 0; j < ncols; j++)
        {
            curr = *breadcrumb_mat_srcentry(mat, i, j);
            top = diag = left = 0;
            if (i) top = *breadcrumb_mat_srcentry(mat, i-1, j);
            if (i && j) diag = *breadcrumb_mat_srcentry(mat, i-1, j-1);
            if (j) left = *breadcrumb_mat_srcentry(mat, i, j-1);

            /*
             * Draw this tableau cell and its connections to its
             * neighbors with non-larger row or column indices.
             */

            /* draw the cell itself, without connections */
            if (curr & CRUMB_CONTENDER) {
                r = 255; g = 0; b = 0;
            } else if (curr & CRUMB_WANT3) {
                r = 0; g = 0; b = 255;
            } else if (curr & CRUMB_WANT2) {
                r = 0; g = 255; b = 0;
            } else {
                r = 0; g = 80; b = 0;
            }
            for (di = 2; di < 5; di++) {
                for (dj = 2; dj < 5; dj++) {
                    buf_set_rgb(buf, j, di, dj, r, g, b);
                }
            }

            /* top connections */
            if ((curr & want_any) && (top & want_any))
            {
                r = 0; b = 0; g = 0;
                if (curr & CRUMB_TOP)
                {
                    r = 0; g = 255; b = 0;
                }
                if ((top & CRUMB_CONTENDER) && (curr & CRUMB_CONTENDER))
                {
                    if (curr & CRUMB_TOP)
                    {
                        r = 255; g = 0; b = 0;
                    }
                }
                buf_set_rgb(buf, j, 0, 3, r, g, b);
                buf_set_rgb(buf, j, 1, 3, r, g, b);
            }

            /* left connections */
            if ((curr & want_any) && (left & want_any))
            {
                r = 0; b = 0; g = 0;
                if (curr & CRUMB_LEFT2)
                {
                    r = 0; g = 255; b = 0;
                }
                if ((left & CRUMB_CONTENDER) && (curr & CRUMB_CONTENDER))
                {
                    if (curr & CRUMB_LEFT)
                    {
                        r = 255; g = 0; b = 0;
                    }
                }
                buf_set_rgb(buf, j, 3, 0, r, g, b);
                buf_set_rgb(buf, j, 3, 1, r, g, b);
            }

            /* diagonal connections */
            if ((curr & want_any) && (diag & want_any))
            {
                r = 0; b = 0; g = 0;
                if (curr & CRUMB_DIAG2)
                {
                    r = 0; g = 255; b = 0;
                }
                if ((diag & CRUMB_CONTENDER) && (curr & CRUMB_CONTENDER))
                {
                    if (curr & CRUMB_DIAG)
                    {
                        r = 255; g = 0; b = 0;
                    }
                }
                buf_set_rgb(buf, j, 0, 0, r, g, b);
                buf_set_rgb(buf, j, 1, 1, r, g, b);
            }
        }

        /* write the row, trimming off the first two rows and cols of pixels */
        for (di = 0; di < 5; di++)
        {
            if (i == 0 && di < 2) continue;
            png_write_row(png_ptr, buf_pixel_row(buf, di) + 2*3);
        }
    }

    png_write_end(png_ptr, NULL);

end:
    if (fout != NULL) fclose(fout);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (info_ptr != NULL) png_destroy_info_struct(png_ptr, &info_ptr);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    buf_clear(buf);

    return code;
}
