/**
 * \file zlib_wrapper.c
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <zlib.h>


int zlib_in_place_compress( void* in, size_t* len )
{
    Byte* compr;
    uLong compr_len;
    int res;

    compr_len = *len + (*len/10 + 1) + 12;

    compr = (Byte*)malloc(compr_len);
    if(compr == NULL) {
        fprintf(stderr, "Cannot allocate memory for zlib compression\n");
        goto error;
    }

    res = compress(compr, &compr_len, in, *len);
    if(res != Z_OK) {
        fprintf(stderr, "Problems compressing with zlib\n");
        goto error;
    }

    if(compr_len >= *len) {
        fprintf(stderr, "zlib didn't compress much\n");
        goto error;
    }

    memcpy(in, compr, compr_len);
    *len = compr_len;

    free(compr);
    return 0;

 error:
    if(compr != NULL) free(compr);
    return -1;
}

int zlib_uncompress( void* out, size_t outlen, const void* in, size_t inlen )
{
    uLong len = outlen;
    int ret = uncompress(out, &len, in, inlen);
    if(ret != Z_OK) {
        fprintf(stderr, "Problems de-compressing with zlib\n");
        return -1;
    }
    
    return 0;
}
