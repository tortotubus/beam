
/* vtkBase64.h
 *
 * A drop‐in “streaming” Base64 encoder for VTK‐style writing.  Instead of
 * calling b64_fwrite(ptr,1,n,stream) multiple times (which emits separate
 * “=”-padded blocks on each call), you can now do:
 *
 *    vtkBase64Writer *ctx = vtkBase64WriterOpen(f);
 *    vtkBase64WriterWrite(ctx, raw1, len1);
 *    vtkBase64WriterWrite(ctx, raw2, len2);
 *    ...
 *    vtkBase64WriterClose(ctx);
 *
 * and you’ll get exactly the same Base64 output as if you had concatenated
 * raw1||raw2||… into one big buffer and called b64_fwrite(...) once.
 */

#include "vtkDataArray.h"

typedef struct vtkBase64Writer_t {
    FILE       *f;
    uint8_t     carry[2];
    int         carry_len;   /* 0..2 */
    char        encode_buf[4];
} vtkBase64Writer;

#if _MPI == 0
/* ----------------------------------------------------------------------
 * 1) Base64 alphabet
 * ---------------------------------------------------------------------- */
static const char _vtk_base64_chars[] =
  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  "abcdefghijklmnopqrstuvwxyz"
  "0123456789+/";

/* ----------------------------------------------------------------------
 * 2) vtkBase64Writer struct
 * 
 *    - `FILE *f`       : the underlying FILE* where we actually fwrite ASCII.
 *    - carry[0..2]     : temporary buffer for up to 2 leftover bytes.
 *    - carry_len       : how many bytes currently in carry (0, 1, or 2).
 *    - encode_buf[4]   : scratch for encoding a 3‐byte group → 4 chars.
 * ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 * 3) Helpers to encode exactly one 3‐byte group → 4 ASCII chars
 * ---------------------------------------------------------------------- */
static void _vtk_base64_encode_triple(const uint8_t *in3, char out4[4])
{
    /* in3[0], in3[1], in3[2] → 24 bits → four 6‐bit chunks */
    uint32_t triple = ((uint32_t)in3[0] << 16)
                    | ((uint32_t)in3[1] <<  8)
                    | ((uint32_t)in3[2]);

    out4[0] = _vtk_base64_chars[(triple >> 18) & 0x3F];
    out4[1] = _vtk_base64_chars[(triple >> 12) & 0x3F];
    out4[2] = _vtk_base64_chars[(triple >>  6) & 0x3F];
    out4[3] = _vtk_base64_chars[(triple      ) & 0x3F];
}

/* ----------------------------------------------------------------------
 * 4) Open a streaming Base64 writer
 * 
 *    Allocates a vtkBase64Writer, zeroes carry buffers, and stores the FILE*.
 *    Returns NULL on allocation failure.
 * ---------------------------------------------------------------------- */
static vtkBase64Writer* vtkBase64WriterOpen(FILE *f)
{
    if (!f) return NULL;
    vtkBase64Writer *ctx = (vtkBase64Writer*)malloc(sizeof(vtkBase64Writer));
    if (!ctx) return NULL;
    ctx->f = f;
    ctx->carry_len = 0;
    return ctx;
}

/* ----------------------------------------------------------------------
 * 5) Write raw bytes into the Base64 stream
 *
 *    - We append the incoming raw bytes to carry[] (0..2 from previous),
 *      so that we can form as many 3‐byte groups as possible.
 *    - Each full 3‐byte group is immediately Base64‐encoded to 4 ASCII
 *      chars (using _vtk_base64_encode_triple) and written via fwrite().
 *    - Any leftover (1 or 2) bytes remain in ctx->carry until the next call.
 *
 *    Returns the number of raw input bytes “consumed” (which will always
 *    equal `size` unless `ctx` is NULL).  If ctx is NULL, returns 0.
 * ---------------------------------------------------------------------- */
static size_t vtkBase64WriterWrite(vtkBase64Writer *ctx,
                                   const void *buf,
                                   size_t size)
{
    if (!ctx || !buf || size == 0) {
        return 0;
    }

    const uint8_t *bytes = (const uint8_t*)buf;
    size_t idx = 0;

    /* If we already have 1 or 2 bytes carried over, try to form a triple */
    if (ctx->carry_len > 0) {
        while (ctx->carry_len < 3 && idx < size) {
            ctx->carry[ctx->carry_len++] = bytes[idx++];
        }
        if (ctx->carry_len == 3) {
            /* We have exactly 3 bytes → encode them */
            _vtk_base64_encode_triple(ctx->carry, ctx->encode_buf);
            fwrite(ctx->encode_buf, 1, 4, ctx->f);
            ctx->carry_len = 0;
        }
    }

    /* Now any remaining full triples in the new buffer: */
    while (idx + 3 <= size) {
        /* encode bytes[idx..idx+2] → 4 chars → write */
        _vtk_base64_encode_triple(&bytes[idx], ctx->encode_buf);
        fwrite(ctx->encode_buf, 1, 4, ctx->f);
        idx += 3;
    }

    /* Any leftover 1 or 2 bytes go into carry[] for next time */
    ctx->carry_len = 0;
    if (idx < size) {
        size_t rem = size - idx;  /* 1 or 2 */
        ctx->carry[0] = bytes[idx];
        if (rem == 2) {
            ctx->carry[1] = bytes[idx + 1];
        }
        ctx->carry_len = rem;
    }

    return size;
}

/* ----------------------------------------------------------------------
 * 6) Close the Base64 stream
 *
 *    If carry_len == 0: do nothing (already aligned to 3).
 *    If carry_len == 1: form a “(a << 16)” triple, encode → two chars + "=="
 *    If carry_len == 2: form a “(a << 16)|(b<<8)” triple, encode → three chars + "="
 *
 *    In each case we emit exactly 4 final ASCII chars.  Then free(ctx).
 *
 *    Returns 0 on success, or -1 if ctx was NULL.
 * ---------------------------------------------------------------------- */
static int vtkBase64WriterClose(vtkBase64Writer *ctx)
{
    if (!ctx) {
        return -1;
    }

    if (ctx->carry_len == 1) {
        /* 1 leftover byte: a = carry[0], triple=(a<<16) → yields 2 chars + “==” */
        uint32_t a = ctx->carry[0];
        uint32_t triple = (a << 16);
        char out4[4];
        out4[0] = _vtk_base64_chars[(triple >> 18) & 0x3F];
        out4[1] = _vtk_base64_chars[(triple >> 12) & 0x3F];
        out4[2] = '=';
        out4[3] = '=';
        fwrite(out4, 1, 4, ctx->f);
    }
    else if (ctx->carry_len == 2) {
        /* 2 leftover bytes: a=carry[0], b=carry[1], triple=(a<<16)|(b<<8) */
        uint32_t a = ctx->carry[0];
        uint32_t b = ctx->carry[1];
        uint32_t triple = (a << 16) | (b << 8);
        char out4[4];
        out4[0] = _vtk_base64_chars[(triple >> 18) & 0x3F];
        out4[1] = _vtk_base64_chars[(triple >> 12) & 0x3F];
        out4[2] = _vtk_base64_chars[(triple >>  6) & 0x3F];
        out4[3] = '=';
        fwrite(out4, 1, 4, ctx->f);
    }
    /* else carry_len == 0 → no final block, everything was a multiple of 3 */

    /* Clean up */
    free(ctx);
    return 0;
}

size_t _vtk_data_array_compute_header_base64(vtkDataArray *vtk_data_array) 
{
    size_t ntuples     = vtk_data_array->number_of_tuples;
    size_t ncomponents = vtk_data_array->number_of_components;
    ncomponents = (ncomponents > 1 ? ncomponents : 1);

    size_t nbytes = ntuples * ncomponents;
    if (vtk_data_array->type == 10) {
        /* bit‐packed: 1 bit per tuple×component, round up to full bytes */
        nbytes = (nbytes + 7) / 8;
    }
    else {
        nbytes *= vtkType_sizeof(vtk_data_array->type);
    }

    return nbytes;
}

#endif