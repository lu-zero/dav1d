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

#include "config.h"

#include "src/mc.h"
#include "src/cpu.h"

#if BITDEPTH == 8 && ARCH_PPC64LE
decl_mc_fn(dav1d_put_8tap_regular_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_regular_smooth_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_regular_sharp_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_smooth_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_smooth_regular_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_smooth_sharp_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_sharp_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_sharp_regular_8bpc_vsx);
decl_mc_fn(dav1d_put_8tap_sharp_smooth_8bpc_vsx);
decl_mc_fn(dav1d_put_bilin_8bpc_vsx);

decl_mct_fn(dav1d_prep_8tap_regular_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_regular_smooth_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_regular_sharp_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_smooth_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_smooth_regular_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_smooth_sharp_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_sharp_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_sharp_regular_8bpc_vsx);
decl_mct_fn(dav1d_prep_8tap_sharp_smooth_8bpc_vsx);
decl_mct_fn(dav1d_prep_bilin_8bpc_vsx);

decl_avg_fn(dav1d_avg_8bpc_vsx);
decl_w_avg_fn(dav1d_w_avg_8bpc_vsx);
decl_mask_fn(dav1d_mask_8bpc_vsx);

decl_warp8x8_fn(dav1d_warp_affine_8x8_8bpc_vsx);
decl_warp8x8t_fn(dav1d_warp_affine_8x8t_8bpc_vsx);
#endif


void bitfn(dav1d_mc_dsp_init_ppc)(Dav1dMCDSPContext *const c) {
#define init_mc_fn(type, name, suffix) \
    c->mc[type] = dav1d_put_##name##_8bpc_##suffix
#define init_mct_fn(type, name, suffix) \
    c->mct[type] = dav1d_prep_##name##_8bpc_##suffix
    const unsigned flags = dav1d_get_cpu_flags();

    if (!(flags & DAV1D_PPC_CPU_FLAG_VSX)) return;

#if BITDEPTH == 8 && ARCH_PPC64LE
    // init_mc_fn (FILTER_2D_8TAP_REGULAR,        8tap_regular,        vsx);
    // init_mc_fn (FILTER_2D_8TAP_REGULAR_SMOOTH, 8tap_regular_smooth, vsx);
    // init_mc_fn (FILTER_2D_8TAP_REGULAR_SHARP,  8tap_regular_sharp,  vsx);
    // init_mc_fn (FILTER_2D_8TAP_SMOOTH_REGULAR, 8tap_smooth_regular, vsx);
    // init_mc_fn (FILTER_2D_8TAP_SMOOTH,         8tap_smooth,         vsx);
    // init_mc_fn (FILTER_2D_8TAP_SMOOTH_SHARP,   8tap_smooth_sharp,   vsx);
    // init_mc_fn (FILTER_2D_8TAP_SHARP_REGULAR,  8tap_sharp_regular,  vsx);
    // init_mc_fn (FILTER_2D_8TAP_SHARP_SMOOTH,   8tap_sharp_smooth,   vsx);
    // init_mc_fn (FILTER_2D_8TAP_SHARP,          8tap_sharp,          vsx);
    // init_mc_fn (FILTER_2D_BILINEAR,            bilin,               vsx);

    // init_mct_fn(FILTER_2D_8TAP_REGULAR,        8tap_regular,        vsx);
    // init_mct_fn(FILTER_2D_8TAP_REGULAR_SMOOTH, 8tap_regular_smooth, vsx);
    // init_mct_fn(FILTER_2D_8TAP_REGULAR_SHARP,  8tap_regular_sharp,  vsx);
    // init_mct_fn(FILTER_2D_8TAP_SMOOTH_REGULAR, 8tap_smooth_regular, vsx);
    // init_mct_fn(FILTER_2D_8TAP_SMOOTH,         8tap_smooth,         vsx);
    // init_mct_fn(FILTER_2D_8TAP_SMOOTH_SHARP,   8tap_smooth_sharp,   vsx);
    // init_mct_fn(FILTER_2D_8TAP_SHARP_REGULAR,  8tap_sharp_regular,  vsx);
    // init_mct_fn(FILTER_2D_8TAP_SHARP_SMOOTH,   8tap_sharp_smooth,   vsx);
    // init_mct_fn(FILTER_2D_8TAP_SHARP,          8tap_sharp,          vsx);
    // init_mct_fn(FILTER_2D_BILINEAR,            bilin,               vsx);

    // c->avg = dav1d_avg_8bpc_vsx;
    // c->w_avg = dav1d_w_avg_8bpc_vsx;
    // c->mask = dav1d_mask_8bpc_vsx;
    // c->warp8x8 = dav1d_warp_affine_8x8_8bpc_vsx;
    // c->warp8x8t = dav1d_warp_affine_8x8t_8bpc_vsx;
#endif
}
