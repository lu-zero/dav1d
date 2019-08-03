/*
 * Copyright © 2019, VideoLAN and dav1d authors
 * Copyright © 2019, Michail Alvanos
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

#include <assert.h>
#include <stdlib.h>

#include "common/intops.h"
#include "src/ppc/types.h"
#include "src/cpu.h"
#include "src/looprestoration.h"

#if BITDEPTH == 8

#define REST_UNIT_STRIDE (390)

static inline i32x4 iclip_vec(i32x4 v, const i32x4 minv, const i32x4 maxv) {
    v = vec_max(minv,v);
    v = vec_min(maxv,v);
    return v;
}

static inline void wiener_filter_h_w256_h64_vsx(uint16_t *hor_ptr,
                                                uint8_t *tmp_ptr,
                                                const int16_t filterh[7])
{
    static const i32x4 zerov = vec_splats( 0 );
    static const i32x4 seven_vec = vec_splats(7);
    static const i32x4 bitdepth_added_vec = vec_splats(1 <<14);
    static const i32x4 round_bits_vec = vec_splats(3);
    static const i32x4 rounding_off_vec = vec_splats(1<<2);
    static const i32x4 clip_limit_v = vec_splats((1 << 13) - 1);

    i32x4 filter_vectors[7] = {
        vec_splats( (int32_t) filterh[0]),
        vec_splats( (int32_t) filterh[1]),
        vec_splats( (int32_t) filterh[2]),
        vec_splats( (int32_t) filterh[3]),
        vec_splats( (int32_t) filterh[4]),
        vec_splats( (int32_t) filterh[5]),
        vec_splats( (int32_t) filterh[6])
    };

    for (int j = 0; j < 64 + 6; j++) {
        for (int i = 0; i < 256; i+=16) {
            i32x4 sum1 = zerov;
            i32x4 sum2 = zerov;
            i32x4 sum3 = zerov;
            i32x4 sum4 = zerov;

            u8x16 tmp_v = vec_vsx_ld(0, &tmp_ptr[i + 3]);

            u16x8 tmp_u16_high = u8h_to_u16(tmp_v);
            u16x8 tmp_u16_low  = u8l_to_u16(tmp_v);

            i32x4 vec_tmp1 = u16h_to_i32(tmp_u16_high);
            i32x4 vec_tmp2 = u16l_to_i32(tmp_u16_high);
            i32x4 vec_tmp3 = u16h_to_i32(tmp_u16_low);
            i32x4 vec_tmp4 = u16l_to_i32(tmp_u16_low);

            for (int k = 0; k < 7; k++) {
                u8x16 ktmp_v = vec_vsx_ld(0, &tmp_ptr[i + k]);

                u16x8 ktmp_u16_high = u8h_to_u16(ktmp_v);
                u16x8 ktmp_u16_low  = u8l_to_u16(ktmp_v);

                i32x4 vec_ktmp1 = u16h_to_i32(ktmp_u16_high);
                i32x4 vec_ktmp2 = u16l_to_i32(ktmp_u16_high);
                i32x4 vec_ktmp3 = u16h_to_i32(ktmp_u16_low);
                i32x4 vec_ktmp4 = u16l_to_i32(ktmp_u16_low);

                sum1 = sum1 + vec_ktmp1 * filter_vectors[k];
                sum2 = sum2 + vec_ktmp2 * filter_vectors[k];
                sum3 = sum3 + vec_ktmp3 * filter_vectors[k];
                sum4 = sum4 + vec_ktmp4 * filter_vectors[k];
            }

            sum1 = sum1 + (vec_tmp1 << seven_vec) + bitdepth_added_vec;
            sum2 = sum2 + (vec_tmp2 << seven_vec) + bitdepth_added_vec;
            sum3 = sum3 + (vec_tmp3 << seven_vec) + bitdepth_added_vec;
            sum4 = sum4 + (vec_tmp4 << seven_vec) + bitdepth_added_vec;

            sum1 = (sum1 + rounding_off_vec) >> round_bits_vec;
            sum2 = (sum2 + rounding_off_vec) >> round_bits_vec;
            sum3 = (sum3 + rounding_off_vec) >> round_bits_vec;
            sum4 = (sum4 + rounding_off_vec) >> round_bits_vec;

            sum1 = iclip_vec(sum1, zerov, clip_limit_v);
            sum2 = iclip_vec(sum2, zerov, clip_limit_v);
            sum3 = iclip_vec(sum3, zerov, clip_limit_v);
            sum4 = iclip_vec(sum4, zerov, clip_limit_v);

            u16x8 sum_short_packed_1 = (u16x8) vec_pack( sum1, sum2 );
            u16x8 sum_short_packed_2 = (u16x8) vec_pack( sum3, sum4 );

            vec_vsx_st(sum_short_packed_1, 0, &hor_ptr[i]);
            vec_vsx_st(sum_short_packed_2, 0, &hor_ptr[i+8]);

        }
        tmp_ptr += REST_UNIT_STRIDE;
        hor_ptr += REST_UNIT_STRIDE;
    }

}

static inline i32x4 iclip_u8_vec(i32x4 v) {
    static const i32x4 zerov = vec_splats(0);
    static const i32x4 maxv = vec_splats(255);

    return iclip_vec(v,zerov,maxv);
}

static inline void wiener_filter_v_w256_h64_vsx(uint8_t *p,
                                                const ptrdiff_t p_stride,
                                                const uint16_t *hor,
                                                const int16_t filterv[7])
{
    static const i32x4 round_bits_vec = vec_splats(11);
    static const i32x4 rounding_off_vec = vec_splats(1<<10);
    static const i32x4 round_offset_vec = vec_splats(1<<18);
    static const i32x4 seven_vec = vec_splats(7);
    static const i32x4 zerov = vec_splats(0);

    i32x4 filter_vectors[7] = {
        vec_splats( (int32_t) filterv[0]),
        vec_splats( (int32_t) filterv[1]),
        vec_splats( (int32_t) filterv[2]),
        vec_splats( (int32_t) filterv[3]),
        vec_splats( (int32_t) filterv[4]),
        vec_splats( (int32_t) filterv[5]),
        vec_splats( (int32_t) filterv[6])
    };

    for (int j = 0; j < 64; j++) {
        for (int i = 0; i < 256; i+=16) {
            u16x8 hor1_v = vec_vsx_ld(0, &hor[(j + 3) * REST_UNIT_STRIDE + i]);
            u16x8 hor2_v = vec_vsx_ld(0, &hor[(j + 3) * REST_UNIT_STRIDE + i + 8]);

            i32x4 hor_v_1 = u16h_to_i32( hor1_v );
            i32x4 hor_v_2 = u16l_to_i32( hor1_v );
            i32x4 hor_v_3 = u16h_to_i32( hor2_v );
            i32x4 hor_v_4 = u16l_to_i32( hor2_v );

            i32x4 sum1 = zerov;
            i32x4 sum2 = zerov;
            i32x4 sum3 = zerov;
            i32x4 sum4 = zerov;

            for (int k = 0; k < 7; k++) {
                u16x8 horf_1 = vec_vsx_ld(0, &hor[(j + k) * REST_UNIT_STRIDE + i]);
                u16x8 horf_2 = vec_vsx_ld(0, &hor[(j + k) * REST_UNIT_STRIDE + i + 8]);

                i32x4 horf_v_1 = u16h_to_i32( horf_1 );
                i32x4 horf_v_2 = u16l_to_i32( horf_1 );
                i32x4 horf_v_3 = u16h_to_i32( horf_2 );
                i32x4 horf_v_4 = u16l_to_i32( horf_2 );

                sum1 = sum1 + horf_v_1 * filter_vectors[k];
                sum2 = sum2 + horf_v_2 * filter_vectors[k];
                sum3 = sum3 + horf_v_3 * filter_vectors[k];
                sum4 = sum4 + horf_v_4 * filter_vectors[k];
            }

            sum1 = (hor_v_1 << seven_vec) - round_offset_vec + sum1;
            sum2 = (hor_v_2 << seven_vec) - round_offset_vec + sum2;
            sum3 = (hor_v_3 << seven_vec) - round_offset_vec + sum3;
            sum4 = (hor_v_4 << seven_vec) - round_offset_vec + sum4;

            sum1 = sum1 + rounding_off_vec;
            sum2 = sum2 + rounding_off_vec;
            sum3 = sum3 + rounding_off_vec;
            sum4 = sum4 + rounding_off_vec;

            sum1 = sum1 >> round_bits_vec;
            sum2 = sum2 >> round_bits_vec;
            sum3 = sum3 >> round_bits_vec;
            sum4 = sum4 >> round_bits_vec;

            sum1 = iclip_u8_vec(sum1);
            sum2 = iclip_u8_vec(sum2);
            sum3 = iclip_u8_vec(sum3);
            sum4 = iclip_u8_vec(sum4);

            u16x8 sum_short_packed_1 = (u16x8) vec_pack( sum1, sum2 );
            u16x8 sum_short_packed_2 = (u16x8) vec_pack( sum3, sum4 );
            u8x16 sum_pixel = vec_pack(sum_short_packed_1, sum_short_packed_2 );

            vec_vsx_st(sum_pixel, 0, &p[j * PXSTRIDE(p_stride) + i]);
        }
    }
}

#define COPY_VEC(d, s, l) \
for (int i=0; i<l; i+=16) \
  vec_vsx_st( vec_vsx_ld( 0, ((uint8_t *)(s)) +i ), 0, ((uint8_t*)(d)) +i) ;

#define COPY_PADDING(d, s, w) \
    COPY_VEC(d, s, w);\
    pixel_copy(d + 256 , s + 256, w - 256);

static inline void padding_vsx(uint8_t *dst, const uint8_t *p,
                               const ptrdiff_t p_stride,
                               const uint8_t (*left)[4],
                               const uint8_t *lpf, const ptrdiff_t lpf_stride,
                               const enum LrEdgeFlags edges)
{
    const int have_left = !!(edges & LR_HAVE_LEFT);
    const int have_right = !!(edges & LR_HAVE_RIGHT);

    int unit_w = 256;
    const int stripe_h = 64;

    // Copy more pixels if we don't have to pad them
    unit_w += 3 * have_left + 3 * have_right;
    uint8_t *dst_l = dst + 3 * !have_left;
    p -= 3 * have_left;
    lpf -= 3 * have_left;


    if (edges & LR_HAVE_TOP) {
        // Copy previous loop filtered rows
        const uint8_t *const above_1 = lpf;
        const uint8_t *const above_2 = above_1 + PXSTRIDE(lpf_stride);

        COPY_PADDING( dst_l, above_1, unit_w);
        COPY_PADDING( dst_l + REST_UNIT_STRIDE, above_1, unit_w);
        COPY_PADDING( dst_l + 2 * REST_UNIT_STRIDE, above_2, unit_w);
    } else {
        // Pad with first row
        COPY_PADDING(dst_l, p, unit_w);
        COPY_PADDING(dst_l + REST_UNIT_STRIDE, p, unit_w);
        COPY_PADDING(dst_l + 2 * REST_UNIT_STRIDE, p, unit_w);
        if (have_left) {
            pixel_copy(dst_l, &left[0][1], 3);
            pixel_copy(dst_l + REST_UNIT_STRIDE, &left[0][1], 3);
            pixel_copy(dst_l + 2 * REST_UNIT_STRIDE, &left[0][1], 3);
        }
    }

    uint8_t *dst_tl = dst_l + 3 * REST_UNIT_STRIDE;
    if (edges & LR_HAVE_BOTTOM) {
        // Copy next loop filtered rows
        const uint8_t *const below_1 = lpf + 6 * PXSTRIDE(lpf_stride);
        const uint8_t *const below_2 = below_1 + PXSTRIDE(lpf_stride);
        COPY_PADDING(dst_tl + stripe_h * REST_UNIT_STRIDE, below_1, unit_w);
        COPY_PADDING(dst_tl + (stripe_h + 1) * REST_UNIT_STRIDE, below_2, unit_w);
        COPY_PADDING(dst_tl + (stripe_h + 2) * REST_UNIT_STRIDE, below_2, unit_w);
    } else {
        // Pad with last row
        const uint8_t *const src = p + (stripe_h - 1) * PXSTRIDE(p_stride);
        COPY_PADDING(dst_tl + stripe_h * REST_UNIT_STRIDE, src, unit_w);
        COPY_PADDING(dst_tl + (stripe_h + 1) * REST_UNIT_STRIDE, src, unit_w);
        COPY_PADDING(dst_tl + (stripe_h + 2) * REST_UNIT_STRIDE, src, unit_w);
        if (have_left) {
            pixel_copy(dst_tl + stripe_h * REST_UNIT_STRIDE, &left[stripe_h - 1][1], 3);
            pixel_copy(dst_tl + (stripe_h + 1) * REST_UNIT_STRIDE, &left[stripe_h - 1][1], 3);
            pixel_copy(dst_tl + (stripe_h + 2) * REST_UNIT_STRIDE, &left[stripe_h - 1][1], 3);
        }
    }

    // Inner UNIT_WxSTRIPE_H
    for (int j = 0; j < stripe_h; j++) {
        COPY_VEC( dst_tl + 3 * have_left, p + 3 * have_left, 240);
        pixel_copy(dst_tl + 3 * have_left + 240, p + 3 * have_left + 240, unit_w - 3 * have_left - 240);
        dst_tl += REST_UNIT_STRIDE;
        p += PXSTRIDE(p_stride);
    }

    if (!have_right) {
        uint8_t *pad = dst_l + unit_w;
        uint8_t *row_last = &dst_l[unit_w - 1];
        // Pad 3x(STRIPE_H+6) with last column
        for (int j = 0; j < stripe_h + 6; j++) {
            pixel_set(pad, *row_last, 3);
            pad += REST_UNIT_STRIDE;
            row_last += REST_UNIT_STRIDE;
        }
    }

    if (!have_left) {
        // Pad 3x(STRIPE_H+6) with first column
        for (int j = 0; j < stripe_h + 6; j++) {
            pixel_set(dst, *dst_l, 3);
            dst += REST_UNIT_STRIDE;
            dst_l += REST_UNIT_STRIDE;
        }
    } else {
        dst += 3 * REST_UNIT_STRIDE;
        for (int j = 0; j < stripe_h; j++) {
            pixel_copy(dst, &left[j][1], 3);
            dst += REST_UNIT_STRIDE;
        }
    }
}

static inline void padding(uint8_t *dst, const uint8_t *p,
                           const ptrdiff_t p_stride, const uint8_t (*left)[4],
                           const uint8_t *lpf, const ptrdiff_t lpf_stride,
                           int unit_w, const int stripe_h,
                           const enum LrEdgeFlags edges)
{
    const int have_left = !!(edges & LR_HAVE_LEFT);
    const int have_right = !!(edges & LR_HAVE_RIGHT);

    // Copy more pixels if we don't have to pad them
    unit_w += 3 * have_left + 3 * have_right;
    uint8_t *dst_l = dst + 3 * !have_left;
    p -= 3 * have_left;
    lpf -= 3 * have_left;

    if (edges & LR_HAVE_TOP) {
        // Copy previous loop filtered rows
        const uint8_t *const above_1 = lpf;
        const uint8_t *const above_2 = above_1 + PXSTRIDE(lpf_stride);
        pixel_copy(dst_l, above_1, unit_w);
        pixel_copy(dst_l + REST_UNIT_STRIDE, above_1, unit_w);
        pixel_copy(dst_l + 2 * REST_UNIT_STRIDE, above_2, unit_w);
    } else {
        // Pad with first row
        pixel_copy(dst_l, p, unit_w);
        pixel_copy(dst_l + REST_UNIT_STRIDE, p, unit_w);
        pixel_copy(dst_l + 2 * REST_UNIT_STRIDE, p, unit_w);
        if (have_left) {
            pixel_copy(dst_l, &left[0][1], 3);
            pixel_copy(dst_l + REST_UNIT_STRIDE, &left[0][1], 3);
            pixel_copy(dst_l + 2 * REST_UNIT_STRIDE, &left[0][1], 3);
        }
    }

    uint8_t *dst_tl = dst_l + 3 * REST_UNIT_STRIDE;
    if (edges & LR_HAVE_BOTTOM) {
        // Copy next loop filtered rows
        const uint8_t *const below_1 = lpf + 6 * PXSTRIDE(lpf_stride);
        const uint8_t *const below_2 = below_1 + PXSTRIDE(lpf_stride);
        pixel_copy(dst_tl + stripe_h * REST_UNIT_STRIDE, below_1, unit_w);
        pixel_copy(dst_tl + (stripe_h + 1) * REST_UNIT_STRIDE, below_2, unit_w);
        pixel_copy(dst_tl + (stripe_h + 2) * REST_UNIT_STRIDE, below_2, unit_w);
    } else {
        // Pad with last row
        const uint8_t *const src = p + (stripe_h - 1) * PXSTRIDE(p_stride);
        pixel_copy(dst_tl + stripe_h * REST_UNIT_STRIDE, src, unit_w);
        pixel_copy(dst_tl + (stripe_h + 1) * REST_UNIT_STRIDE, src, unit_w);
        pixel_copy(dst_tl + (stripe_h + 2) * REST_UNIT_STRIDE, src, unit_w);
        if (have_left) {
            pixel_copy(dst_tl + stripe_h * REST_UNIT_STRIDE, &left[stripe_h - 1][1], 3);
            pixel_copy(dst_tl + (stripe_h + 1) * REST_UNIT_STRIDE, &left[stripe_h - 1][1], 3);
            pixel_copy(dst_tl + (stripe_h + 2) * REST_UNIT_STRIDE, &left[stripe_h - 1][1], 3);
        }
    }

    // Inner UNIT_WxSTRIPE_H
    for (int j = 0; j < stripe_h; j++) {
        pixel_copy(dst_tl + 3 * have_left, p + 3 * have_left, unit_w - 3 * have_left);
        dst_tl += REST_UNIT_STRIDE;
        p += PXSTRIDE(p_stride);
    }

    if (!have_right) {
        uint8_t *pad = dst_l + unit_w;
        uint8_t *row_last = &dst_l[unit_w - 1];
        // Pad 3x(STRIPE_H+6) with last column
        for (int j = 0; j < stripe_h + 6; j++) {
            pixel_set(pad, *row_last, 3);
            pad += REST_UNIT_STRIDE;
            row_last += REST_UNIT_STRIDE;
        }
    }

    if (!have_left) {
        // Pad 3x(STRIPE_H+6) with first column
        for (int j = 0; j < stripe_h + 6; j++) {
            pixel_set(dst, *dst_l, 3);
            dst += REST_UNIT_STRIDE;
            dst_l += REST_UNIT_STRIDE;
        }
    } else {
        dst += 3 * REST_UNIT_STRIDE;
        for (int j = 0; j < stripe_h; j++) {
            pixel_copy(dst, &left[j][1], 3);
            dst += REST_UNIT_STRIDE;
        }
    }
}

static inline void wiener_filter_h(uint16_t *hor_ptr, uint8_t *tmp_ptr,
                                   const int16_t filterh[7], const int w,
                                   const int h)
{
    const int bitdepth = 8;
    const int round_bits_h = 3 + (bitdepth == 12) * 2;
    const int rounding_off_h = 1 << (round_bits_h - 1);
    const int clip_limit = 1 << (bitdepth + 1 + 7 - round_bits_h);
    for (int j = 0; j < h + 6; j++) {
        for (int i = 0; i < w; i++) {
            int sum = (tmp_ptr[i + 3] << 7) + (1 << (bitdepth + 6));

            for (int k = 0; k < 7; k++) {
                sum += tmp_ptr[i + k] * filterh[k];
            }

            hor_ptr[i] =
                iclip((sum + rounding_off_h) >> round_bits_h, 0, clip_limit - 1);
        }
        tmp_ptr += REST_UNIT_STRIDE;
        hor_ptr += REST_UNIT_STRIDE;
    }
}

static inline void wiener_filter_v(uint8_t *p, const ptrdiff_t p_stride,
                                   const uint16_t *hor,
                                   const int16_t filterv[7],
                                   const int w, const int h)
{
    const int bitdepth = 8;
    const int round_bits_v = 11 - (bitdepth == 12) * 2;
    const int rounding_off_v = 1 << (round_bits_v - 1);
    const int round_offset = 1 << (bitdepth + (round_bits_v - 1));
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            int sum = (hor[(j + 3) * REST_UNIT_STRIDE + i] << 7) - round_offset;

            for (int k = 0; k < 7; k++) {
                sum += hor[(j + k) * REST_UNIT_STRIDE + i] * filterv[k];
            }

            p[j * PXSTRIDE(p_stride) + i] =
                iclip_pixel((sum + rounding_off_v) >> round_bits_v);
        }
    }
}

// FIXME Could split into luma and chroma specific functions,
// (since first and last tops are always 0 for chroma)
// FIXME Could implement a version that requires less temporary memory
// (should be possible to implement with only 6 rows of temp storage)
static void wiener_filter_vsx(uint8_t *p, const ptrdiff_t p_stride,
                              const uint8_t (*const left)[4],
                              const uint8_t *lpf,
                              const ptrdiff_t lpf_stride,
                              const int w, const int h,
                              const int16_t filterh[7],
                              const int16_t filterv[7],
                              const enum LrEdgeFlags edges HIGHBD_DECL_SUFFIX)
{
    // Wiener filtering is applied to a maximum stripe height of 64 + 3 pixels
    // of padding above and below
    ALIGN(uint8_t tmp[70 /*(64 + 3 + 3)*/ * REST_UNIT_STRIDE], 16);
    ALIGN(uint16_t hor[70 /*(64 + 3 + 3)*/ * REST_UNIT_STRIDE], 16);

    if (w == 256 && h == 64){
        padding_vsx(tmp, p, p_stride, left, lpf, lpf_stride, edges);
        wiener_filter_h_w256_h64_vsx(hor, tmp, filterh);
        wiener_filter_v_w256_h64_vsx(p, p_stride, hor, filterv);
    } else{

        padding(tmp, p, p_stride, left, lpf, lpf_stride, w, h, edges);
        wiener_filter_h(hor, tmp, filterh, w, h);
        wiener_filter_v(p, p_stride, hor, filterv, w, h);
    }
}
#endif

COLD void bitfn(dav1d_loop_restoration_dsp_init_ppc)
    (Dav1dLoopRestorationDSPContext *const c)
{
    const unsigned flags = dav1d_get_cpu_flags();

    if (!(flags & DAV1D_PPC_CPU_FLAG_VSX)) return;

#if BITDEPTH == 8
    c->wiener = wiener_filter_vsx;
#endif
}


