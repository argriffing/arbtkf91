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

#include "dp.h"
#include "vis.h"

/* RGBA */
#define PIXEL_WIDTH 4

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
static void buf_set_rgba(buf_t buf, slong j, slong di, slong dj,
        png_byte r, png_byte g, png_byte b, png_byte a);

void
buf_init(buf_t buf, slong nrows, slong ncols)
{
    buf->sz_data = 5 * ncols * 5 * PIXEL_WIDTH * sizeof(png_byte);
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
    return buf->data + di * (buf->c * 5 * PIXEL_WIDTH);
}

png_bytep
buf_pixel(buf_t buf, slong j, slong di, slong dj)
{
    /* the buffer has rgb values for each 5 pixel x 5 pixel cell
     * in one row of the tableau.
     * j is the column of the dp cell
     * di and dj index the row and column of the pixel within the dp cell */
    return buf_pixel_row(buf, di) + j * (5 * PIXEL_WIDTH) + dj * PIXEL_WIDTH;
}

void
buf_set_rgba(buf_t buf, slong j, slong di, slong dj,
        png_byte r, png_byte g, png_byte b, png_byte a)
{
    png_bytep pixel;
    pixel = buf_pixel(buf, j, di, dj);
    pixel[0] = r;
    pixel[1] = g;
    pixel[2] = b;
    pixel[3] = a;
}

void
buf_zero(buf_t buf)
{
    memset(buf->data, 0, buf->sz_data);
}

static __inline__ void
_firebrick(png_byte *r, png_byte *g, png_byte *b, png_byte *a)
{
    *r = 0xB2; *g = 0x22; *b = 0x22; *a = 0xFF;
}

static __inline__ void
_myblue(png_byte *r, png_byte *g, png_byte *b, png_byte *a)
{
    *r = 0x1E; *g = 0x36; *b = 0x91; *a = 0xFF;
}


int write_tableau_image(const char * filename,
        dp_mat_t mat, const char * title)
{
    FILE *fout;
    png_structp png_ptr;
    png_infop info_ptr;
    int code, width, height;
    png_byte r, g, b, a;
    dp_t curr;

    fout = NULL;
    png_ptr = NULL;
    info_ptr = NULL;
    code = 0;

    slong nrows, ncols;

    nrows = dp_mat_nrows(mat);
    ncols = dp_mat_ncols(mat);

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
            8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != NULL)
    {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = "Title";
        title_text.text = (char *) title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    /* write image data one row at a time */
    int i, j, di, dj;
    buf_t buf;
    buf_init(buf, nrows, ncols);
    for (i = 0; i < nrows; i++)
    {
        /* reset all entries of the pixel buffer to zero */
        buf_zero(buf);

        for (j = 0; j < ncols; j++)
        {
            curr = *dp_mat_entry(mat, i, j);

            /*
             * Draw this tableau cell and its connections to its
             * neighbors with non-larger row or column indices.
             *
             * Note that it's not important to be careful about
             * drawing connections from tableau cells on the boundary,
             * because those boundary pixels will be written into the png.
             */

            /* draw the cell itself, without connections */
            r = 0; g = 0; b = 0; a = 0;
            if (curr & DP_TRACE) {
                _myblue(&r, &g, &b, &a);
            } else if (curr & DP_MAX3) {
                _firebrick(&r, &g, &b, &a);
            } else if (curr & DP_MAX2) {
                r = 0; g = 255; b = 0; a = 255;
            } else {
                r = 245; g = 245; b = 245; a = 255;
            }
            for (di = 2; di < 5; di++) {
                for (dj = 2; dj < 5; dj++) {
                    buf_set_rgba(buf, j, di, dj, r, g, b, a);
                }
            }

            /* top connection */
            if (dp_m0_is_interesting(curr))
            {
                r = 0; g = 255; b = 0; a=255;
                if ((curr & DP_TRACE) && (curr & DP_MAX3_M0))
                {
                    _myblue(&r, &g, &b, &a);
                }
                buf_set_rgba(buf, j, 0, 3, r, g, b, a);
                buf_set_rgba(buf, j, 1, 3, r, g, b, a);
            }

            /* diagonal connection */
            if (dp_m1_is_interesting(curr))
            {
                r = 0; g = 255; b = 0; a=255;
                if ((curr & DP_TRACE) && (curr & DP_MAX3_M1))
                {
                    _myblue(&r, &g, &b, &a);
                }
                buf_set_rgba(buf, j, 0, 0, r, g, b, a);
                buf_set_rgba(buf, j, 1, 1, r, g, b, a);
            }

            /* left connection */
            if (dp_m2_is_interesting(curr))
            {
                r = 0; g = 255; b = 0; a=255;
                if ((curr & DP_TRACE) && (curr & DP_MAX3_M2))
                {
                    _myblue(&r, &g, &b, &a);
                }
                buf_set_rgba(buf, j, 3, 0, r, g, b, a);
                buf_set_rgba(buf, j, 3, 1, r, g, b, a);
            }

        }

        /* write the row, trimming off the first two rows and cols of pixels */
        for (di = 0; di < 5; di++)
        {
            if (i == 0 && di < 2) continue;
            png_write_row(png_ptr, buf_pixel_row(buf, di) + 2*PIXEL_WIDTH);
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


int write_simple_tableau_image(const char * filename,
        dp_mat_t mat, const char * title)
{
    FILE *fout;
    png_structp png_ptr;
    png_infop info_ptr;
    int code, width, height;
    png_byte r, g, b, a;
    dp_t curr;
    png_bytep pixel_row;

    fout = NULL;
    png_ptr = NULL;
    info_ptr = NULL;
    code = 0;

    slong nrows, ncols;

    nrows = dp_mat_nrows(mat);
    ncols = dp_mat_ncols(mat);

    width = (int) ncols;
    height = (int) nrows;

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
            8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != NULL)
    {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = "Title";
        title_text.text = (char *) title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    /* write image data one row at a time */
    size_t sz_pixel_row = width * PIXEL_WIDTH * sizeof(png_byte);
    pixel_row = malloc(sz_pixel_row);
    int i, j;
    for (i = 0; i < nrows; i++)
    {
        /* reset all entries of the pixel buffer to zero */
        memset(pixel_row, 0, sz_pixel_row);

        for (j = 0; j < ncols; j++)
        {
            curr = *dp_mat_entry(mat, i, j);

            /* draw the cell itself, without connections */
            r = 0; g = 0; b = 0; a = 0;
            if (curr & DP_TRACE) {
                _myblue(&r, &g, &b, &a);
            } else if (curr & DP_MAX3) {
                _firebrick(&r, &g, &b, &a);
            } else if (curr & DP_MAX2) {
                r = 0; g = 255; b = 0; a = 255;
            }
            pixel_row[PIXEL_WIDTH * j + 0] = r;
            pixel_row[PIXEL_WIDTH * j + 1] = g;
            pixel_row[PIXEL_WIDTH * j + 2] = b;
            pixel_row[PIXEL_WIDTH * j + 3] = a;
        }

        /* write the row */
        png_write_row(png_ptr, pixel_row);
    }

    png_write_end(png_ptr, NULL);

end:
    if (fout != NULL) fclose(fout);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (info_ptr != NULL) png_destroy_info_struct(png_ptr, &info_ptr);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

    free(pixel_row);

    return code;
}
