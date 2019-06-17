/*
 * Copyright © 2019, VideoLAN and dav1d authors
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <altivec.h>

/***********************************************************************
 * Vector types
 **********************************************************************/
#define vec_u8_t  vector unsigned char
#define vec_s8_t  vector signed char
#define vec_u16_t vector unsigned short
#define vec_s16_t vector signed short
#define vec_u32_t vector unsigned int
#define vec_s32_t vector signed int

#ifdef __VSX__
#define vec_u64_t vector unsigned long long
#define vec_s64_t vector signed long long
#endif

/***********************************************************************
 * Null vector
 **********************************************************************/
#define LOAD_ZERO const vec_u8_t zerov = vec_splat_u8( 0 )

#define zero_u8v  (vec_u8_t)  zerov
#define zero_s8v  (vec_s8_t)  zerov
#define zero_u16v (vec_u16_t) zerov
#define zero_s16v (vec_s16_t) zerov
#define zero_u32v (vec_u32_t) zerov
#define zero_s32v (vec_s32_t) zerov

/***********************************************************************
 * 8 <-> 16 bits conversions
 **********************************************************************/
#ifdef WORDS_BIGENDIAN
#define vec_u8_to_u16_h(v) (vec_u16_t) vec_mergeh( zero_u8v, (vec_u8_t) v )
#define vec_u8_to_u16_l(v) (vec_u16_t) vec_mergel( zero_u8v, (vec_u8_t) v )
#define vec_u8_to_s16_h(v) (vec_s16_t) vec_mergeh( zero_u8v, (vec_u8_t) v )
#define vec_u8_to_s16_l(v) (vec_s16_t) vec_mergel( zero_u8v, (vec_u8_t) v )
#else
#define vec_u8_to_u16_h(v) (vec_u16_t) vec_mergeh( (vec_u8_t) v, zero_u8v )
#define vec_u8_to_u16_l(v) (vec_u16_t) vec_mergel( (vec_u8_t) v, zero_u8v )
#define vec_u8_to_s16_h(v) (vec_s16_t) vec_mergeh( (vec_u8_t) v, zero_u8v )
#define vec_u8_to_s16_l(v) (vec_s16_t) vec_mergel( (vec_u8_t) v, zero_u8v )
#endif

#define vec_u8_to_u16(v) vec_u8_to_u16_h(v)
#define vec_u8_to_s16(v) vec_u8_to_s16_h(v)

#define vec_u16_to_u8(v) vec_pack( v, zero_u16v )
#define vec_s16_to_u8(v) vec_packsu( v, zero_s16v )


/***********************************************************************
 * 16 <-> 32 bits conversions
 **********************************************************************/
#ifdef WORDS_BIGENDIAN
#define vec_u16_to_u32_h(v) (vec_u32_t) vec_mergeh( zero_u16v, (vec_u16_t) v )
#define vec_u16_to_u32_l(v) (vec_u32_t) vec_mergel( zero_u16v, (vec_u16_t) v )
#define vec_u16_to_s32_h(v) (vec_s32_t) vec_mergeh( zero_u16v, (vec_u16_t) v )
#define vec_u16_to_s32_l(v) (vec_s32_t) vec_mergel( zero_u16v, (vec_u16_t) v )
#else
#define vec_u16_to_u32_h(v) (vec_u32_t) vec_mergeh( (vec_u16_t) v, zero_u16v )
#define vec_u16_to_u32_l(v) (vec_u32_t) vec_mergel( (vec_u16_t) v, zero_u16v )
#define vec_u16_to_s32_h(v) (vec_s32_t) vec_mergeh( (vec_u16_t) v, zero_u16v )
#define vec_u16_to_s32_l(v) (vec_s32_t) vec_mergel( (vec_u16_t) v, zero_u16v )
#endif

#define vec_u16_to_u32(v) vec_u16_to_u32_h(v)
#define vec_u16_to_s32(v) vec_u16_to_s32_h(v)

#define vec_u32_to_u16(v) vec_pack( v, zero_u32v )
#define vec_s32_to_u16(v) vec_packsu( v, zero_s32v )


/***********************************************************************
 * Vector Loads and Stores
 ***********************************************************************/

#ifndef __VSX__

#undef vec_vsx_ld
#define vec_vsx_ld(off, src) \
    vec_perm(vec_ld(off, src), vec_ld(off + 15, src), vec_lvsl(off, src))

#undef vec_vsx_st
#define vec_vsx_st(v, off, dst)                            \
    do {                                                   \
        uint8_t *_dst = (uint8_t*)(dst);                   \
        vec_u8_t _v = (vec_u8_t)(v);                       \
        vec_u8_t _a = vec_ld(off, _dst);                   \
        vec_u8_t _b = vec_ld(off + 15, _dst);              \
        vec_u8_t _e = vec_perm(_b, _a, vec_lvsl(0, _dst)); \
        vec_u8_t _m = vec_lvsr(0, _dst);                   \
                                                           \
        vec_st(vec_perm(_v, _e, _m), off + 15, _dst);      \
        vec_st(vec_perm(_e, _v, _m), off, _dst);           \
    } while( 0 )

#endif


/***********************************************************************
 * Vector Loads and Stores
 ***********************************************************************/

#ifndef __POWER9_VECTOR__
#define tiny_copy( d, s, l ) memcpy( d, s, l )
#else
#define tiny_copy( d, s, l ) vec_xst_len( vec_vsx_ld( 0, s ), d, l )
#endif

#define COPY16(src,dst) vec_vsx_st(vec_vsx_ld(0, src), 0, dst);


/***********************************************************************
 * Data aligment
 ***********************************************************************/

#define DECLARE_ALIGNED( var, n ) var __attribute__((aligned(n)))

#define ALIGNED_4( var )  DECLARE_ALIGNED( var, 4 )
#define ALIGNED_8( var )  DECLARE_ALIGNED( var, 8 )
#define ALIGNED_16( var ) DECLARE_ALIGNED( var, 16 )



