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

#if BITDEPTH == 8 && ARCH_PPC64LE

void sgr_filter_vsx(pixel *const dst, const ptrdiff_t dst_stride,
                    const pixel (*const left)[4],
                    const pixel *lpf, const ptrdiff_t lpf_stride,
                    const int w, const int h, const int sgr_idx,
                    const int16_t sgr_wt[7], const enum LrEdgeFlags edges);

void wiener_filter_vsx(pixel *const dst, const ptrdiff_t dst_stride,
                       const pixel (*const left)[4],
                       const pixel *lpf, const ptrdiff_t lpf_stride,
                       const int w, const int h, const int16_t fh[7],
                       const int16_t fv[7], const enum LrEdgeFlags edges);

#endif

COLD void bitfn(dav1d_loop_restoration_dsp_init_ppc)(Dav1dLoopRestorationDSPContext *const c) {
    const unsigned flags = dav1d_get_cpu_flags();

    if (!(flags & DAV1D_PPC_CPU_FLAG_VSX)) return;

#if BITDEPTH == 8 && ARCH_PPC64LE
    // c->selfguided = sgr_filter_vsx;
    // c->wiener = wiener_filter_vsx;
#endif
}

