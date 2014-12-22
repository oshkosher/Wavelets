/**
 * \file zlib_wrapper.h
 * 
 */
#ifndef __ZLIB_WRAPPER_H__
#define __ZLIB_WRAPPER_H__

#include <zlib.h>

int zlib_in_place_compress( void* in, size_t* len );

int zlib_uncompress( void* out, size_t outlen, const void* in, size_t inlen );


#endif /* __ZLIB_WRAPPER_H__ */
