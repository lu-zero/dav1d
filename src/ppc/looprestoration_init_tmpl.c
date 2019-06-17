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

#include "src/cpu.h"
#include "src/looprestoration.h"
#include "common/intops.h"
#include "ppc_common.h"

// 256 * 1.5 + 3 + 3 = 390
#define REST_UNIT_STRIDE (390)

#define INLINE inline

#if BITDEPTH == 8 && ARCH_PPC64LE

void sgr_filter_vsx(pixel *const dst, const ptrdiff_t dst_stride,
                    const pixel (*const left)[4],
                    const pixel *lpf, const ptrdiff_t lpf_stride,
                    const int w, const int h, const int sgr_idx,
                    const int16_t sgr_wt[7], const enum LrEdgeFlags edges);


static INLINE vec_s32_t iclip_vec(vec_s32_t v, const int min, const int max) {
    vec_s32_t minv = vec_splats( (int) min );
    vec_s32_t maxv = vec_splats( (int) max );
    v = vec_max(minv,v);
    v = vec_min(maxv,v);
    return v;
}


static INLINE  void wiener_filter_h_w256_h64_vsx(
        uint16_t *hor_ptr,  uint8_t *tmp_ptr,
        const int16_t filterh[7],
        const int bitdepth
        ){

    const int round_bits_h = 3 + (bitdepth == 12) * 2;
    const int rounding_off_h = 1 << (round_bits_h - 1);
    const int clip_limit = 1 << (bitdepth + 1 + 7 - round_bits_h);

    LOAD_ZERO;

    vec_s32_t round_bits_vec = vec_splats( (signed int) round_bits_h);
    vec_s32_t rounding_off_vec = vec_splats(rounding_off_h);
    vec_s32_t seven_vec = vec_splats(7);
    vec_s32_t bitdepth_added_vec = vec_splats(1 << (bitdepth + 6));

    vec_s32_t filter_vectors[7] = {
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

            vec_s32_t sum1 = zero_s32v;
            vec_s32_t sum2 = zero_s32v;
            vec_s32_t sum3 = zero_s32v;
            vec_s32_t sum4 = zero_s32v;

            vec_u8_t tmp_v = vec_vsx_ld(0, &tmp_ptr[i + 3]);

            vec_u16_t tmp_u16_high = vec_u8_to_u16_h(tmp_v);
            vec_u16_t tmp_u16_low  = vec_u8_to_u16_l(tmp_v);

            vec_s32_t vec_tmp1 = vec_u16_to_s32_h(tmp_u16_high);
            vec_s32_t vec_tmp2 = vec_u16_to_s32_l(tmp_u16_high);
            vec_s32_t vec_tmp3 = vec_u16_to_s32_h(tmp_u16_low);
            vec_s32_t vec_tmp4 = vec_u16_to_s32_l(tmp_u16_low);

            for (int k = 0; k < 7; k++) {
                vec_u8_t ktmp_v = vec_vsx_ld(0, &tmp_ptr[i + k]);

                vec_u16_t ktmp_u16_high = vec_u8_to_u16_h(ktmp_v);
                vec_u16_t ktmp_u16_low  = vec_u8_to_u16_l(ktmp_v);

                vec_s32_t vec_tmp1 = vec_u16_to_s32_h(ktmp_u16_high);
                vec_s32_t vec_tmp2 = vec_u16_to_s32_l(ktmp_u16_high);
                vec_s32_t vec_tmp3 = vec_u16_to_s32_h(ktmp_u16_low);
                vec_s32_t vec_tmp4 = vec_u16_to_s32_l(ktmp_u16_low);

                sum1 = sum1 + vec_tmp1 * filter_vectors[k];
                sum2 = sum2 + vec_tmp2 * filter_vectors[k];
                sum3 = sum3 + vec_tmp3 * filter_vectors[k];
                sum4 = sum4 + vec_tmp4 * filter_vectors[k];
            }

            sum1 = sum1 + ( vec_tmp1 << seven_vec) + bitdepth_added_vec;
            sum2 = sum2 + ( vec_tmp2 << seven_vec) + bitdepth_added_vec;
            sum3 = sum3 + ( vec_tmp3 << seven_vec) + bitdepth_added_vec;
            sum4 = sum4 + ( vec_tmp4 << seven_vec) + bitdepth_added_vec;

            sum1 = (sum1 + rounding_off_vec) >> round_bits_vec;
            sum2 = (sum2 + rounding_off_vec) >> round_bits_vec;
            sum3 = (sum3 + rounding_off_vec) >> round_bits_vec;
            sum4 = (sum4 + rounding_off_vec) >> round_bits_vec;

            sum1 = iclip_vec(sum1, 0, clip_limit -1);
            sum2 = iclip_vec(sum2, 0, clip_limit -1);
            sum3 = iclip_vec(sum3, 0, clip_limit -1);
            sum4 = iclip_vec(sum4, 0, clip_limit -1);

            vec_u16_t sum_short_packed_1 = (vec_u16_t) vec_pack( sum1, sum2 );
            vec_u16_t sum_short_packed_2 = (vec_u16_t) vec_pack( sum3, sum4 );

            vec_vsx_st(sum_short_packed_1, 0, &hor_ptr[i]);
            vec_vsx_st(sum_short_packed_2, 0, &hor_ptr[i+8]);


        }
        tmp_ptr += REST_UNIT_STRIDE;
        hor_ptr += REST_UNIT_STRIDE;
    }

}

static inline vec_s32_t iclip_u8_vec(vec_s32_t v) {
    return iclip_vec(v,0,255);
}

static  INLINE  void wiener_filter_v_w256_h64_vsx(
        uint8_t *p,
        const ptrdiff_t p_stride, const uint16_t *hor,
        const int16_t filterv[7],
        const int bitdepth
        ){
    const int round_bits_v = 11 - (bitdepth == 12) * 2;
    const int rounding_off_v = 1 << (round_bits_v - 1);
    const int round_offset = 1 << (bitdepth + (round_bits_v - 1));

    vec_s32_t round_bits_vec = vec_splats( (signed int) round_bits_v);
    vec_s32_t rounding_off_vec = vec_splats(rounding_off_v);
    vec_s32_t round_offset_vec = vec_splats(round_offset);
    vec_s32_t seven_vec = vec_splats(7);

    LOAD_ZERO;
    vec_s32_t filter_vectors[7] = {
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

            vec_s32_t sum1 = zero_s32v;
            vec_s32_t sum2 = zero_s32v;
            vec_s32_t sum3 = zero_s32v;
            vec_s32_t sum4 = zero_s32v;

            vec_u16_t hor1_v = vec_vsx_ld(0, &hor[(j + 3) * REST_UNIT_STRIDE + i]);
            vec_u16_t hor2_v = vec_vsx_ld(0, &hor[(j + 3) * REST_UNIT_STRIDE + i + 8]);

            vec_s32_t hor_high_1 = vec_u16_to_s32_h(hor1_v);
            vec_s32_t hor_low_1  = vec_u16_to_s32_l(hor1_v);

            vec_s32_t hor_low_2  = vec_u16_to_s32_l(hor2_v);
            vec_s32_t hor_high_2 = vec_u16_to_s32_h(hor2_v);

            for (int k = 0; k < 7; k++) {

                vec_u16_t hor_v_1 = vec_vsx_ld(0, &hor[(j + k) * REST_UNIT_STRIDE + i]);
                vec_s32_t hor_high_1 = vec_u16_to_s32_h(hor_v_1);
                vec_s32_t hor_low_1  = vec_u16_to_s32_l(hor_v_1);

                vec_u16_t hor_v_2 = vec_vsx_ld(0, &hor[(j + k) * REST_UNIT_STRIDE + i + 8]);
                vec_s32_t hor_high_2 = vec_u16_to_s32_h(hor_v_2);
                vec_s32_t hor_low_2  = vec_u16_to_s32_l(hor_v_2);

                sum1 = sum1 + hor_high_1 * filter_vectors[k];
                sum2 = sum2 + hor_low_1  * filter_vectors[k];
                sum3 = sum3 + hor_high_2 * filter_vectors[k];
                sum4 = sum4 + hor_low_2  * filter_vectors[k];
            }

            sum1 = sum1 + (hor_high_1 << seven_vec) - round_offset_vec;
            sum2 = sum2 + (hor_low_1  << seven_vec) - round_offset_vec;
            sum3 = sum3 + (hor_high_2 << seven_vec) - round_offset_vec;
            sum4 = sum4 + (hor_low_2  << seven_vec) - round_offset_vec;

            sum1 = (sum1 + rounding_off_vec) >> round_bits_vec;
            sum2 = (sum2 + rounding_off_vec) >> round_bits_vec;
            sum3 = (sum3 + rounding_off_vec) >> round_bits_vec;
            sum4 = (sum4 + rounding_off_vec) >> round_bits_vec;

            sum1 = iclip_u8_vec(sum1);
            sum2 = iclip_u8_vec(sum2);
            sum3 = iclip_u8_vec(sum3);
            sum4 = iclip_u8_vec(sum4);

            vec_u16_t sum_short_packed_1 = (vec_u16_t) vec_pack( sum1, sum2 );
            vec_u16_t sum_short_packed_2 = (vec_u16_t) vec_pack( sum3, sum4 );
            vec_u8_t sum_pixel = vec_pack(sum_short_packed_1, sum_short_packed_2 );

            vec_vsx_st(sum_pixel, 0, &p[j * PXSTRIDE(p_stride) + i]);
        }
    }
}



/// Generic

// TODO Reuse p when no padding is needed (add and remove lpf pixels in p)
// TODO Chroma only requires 2 rows of padding.
static   INLINE  void padding_vsx(uint8_t *dst, const uint8_t *p, const ptrdiff_t p_stride,
                    const uint8_t (*left)[4],
                    const uint8_t *lpf, const ptrdiff_t lpf_stride,
                    int unit_w, const int stripe_h, const enum LrEdgeFlags edges)
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


static  INLINE  void wiener_filter_h_vsx(
        uint16_t *hor_ptr,  uint8_t *tmp_ptr,
        const int16_t filterh[7],
        const int w, const int h,
        const int bitdepth
        ){

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



static INLINE void wiener_filter_v_vsx(
        uint8_t *p,
        const ptrdiff_t p_stride, const uint16_t *hor,
        const int16_t filterv[7],
        const int w, const int h,
        const int bitdepth
        ){
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
static INLINE void wiener_filter_vsx(uint8_t *p, const ptrdiff_t p_stride,
                     const uint8_t (*const left)[4],
                     const uint8_t *lpf, const ptrdiff_t lpf_stride,
                     const int w, const int h,
                     const int16_t filterh[7], const int16_t filterv[7],
                     const enum LrEdgeFlags edges HIGHBD_DECL_SUFFIX)
{
    // Wiener filtering is applied to a maximum stripe height of 64 + 3 pixels
    // of padding above and below
    ALIGNED_16(uint8_t tmp[70 /*(64 + 3 + 3)*/ * REST_UNIT_STRIDE]);

    // Values stored between horizontal and vertical filtering don't
    // fit in a uint8_t.
    ALIGNED_16(uint16_t hor[70 /*(64 + 3 + 3)*/ * REST_UNIT_STRIDE]);

    const int bitdepth = bitdepth_from_max(bitdepth_max);

    padding_vsx(tmp, p, p_stride, left, lpf, lpf_stride, w, h, edges);

    if (w == 256 && h == 64){
        wiener_filter_h_w256_h64_vsx(hor, tmp, filterh, bitdepth);
        wiener_filter_v_w256_h64_vsx(p, p_stride, hor, filterv, bitdepth  );
    } else{
        wiener_filter_h_vsx(hor, tmp, filterh, w, h, bitdepth);
        wiener_filter_v_vsx(p, p_stride, hor, filterv, w, h, bitdepth  );
    }

}





#endif

COLD void bitfn(dav1d_loop_restoration_dsp_init_ppc)(Dav1dLoopRestorationDSPContext *const c) {
    const unsigned flags = dav1d_get_cpu_flags();

    if (!(flags & DAV1D_PPC_CPU_FLAG_VSX)) return;

#if BITDEPTH == 8 && ARCH_PPC64LE
    // c->selfguided = sgr_filter_vsx;
    c->wiener = wiener_filter_vsx;
#endif
}

