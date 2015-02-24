/*
 * Copyright (c) 2003-2010 University of Florida
 * Copyright (c) 2013-2015 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * The GNU Lesser General Public License is included in this distribution
 * in the file COPYING.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "boys.h"
#include "erd.h"
#include "erdutil.h"

#ifdef __x86_64__
#include <x86intrin.h>
#elif defined(__MIC__)
#include <immintrin.h>
#endif

ERD_OFFLOAD static inline double pow3o4(double x) {
    return __builtin_sqrt(x * __builtin_sqrt(x));
}

ERD_OFFLOAD static inline double square(double x) {
    return x * x;
}

ERD_OFFLOAD static const double TOL = 1.0e-14;

ERD_OFFLOAD static const double c0 = 0x1.0B1A240FD5AF4p-8;
ERD_OFFLOAD static const double c1 = 0x1.352866F31ED93p+0;
ERD_OFFLOAD static const double c2 = -0x1.567450B98A180p-1;
ERD_OFFLOAD static const double c3 = 0x1.9721DD0ADF393p-3;
ERD_OFFLOAD static const double c4 = -0x1.EF4155EFB6D81p-6;
ERD_OFFLOAD static const double c5 = 0x1.DFF4D0D064A26p-10;
ERD_OFFLOAD static const double x0 = 0x1.628C5E7D820BFp-5;

#ifdef __AVX__
    __attribute__((aligned(16))) static const uint8_t xmm_pack_table[16*16] = {
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7,    8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80,
          12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80,
           8,    9,   10,   11,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    8,    9,   10,   11,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15
    };

    static const __m256d ymm_c0 = {  0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8  };
    static const __m256d ymm_c1 = {  0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0  };
    static const __m256d ymm_c2 = { -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1  };
    static const __m256d ymm_c3 = {  0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3  };
    static const __m256d ymm_c4 = { -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6  };
    static const __m256d ymm_c5 = {  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10 };
    static const __m256d ymm_x0 = {  0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5  };
    static const __m256d ymm_one = { 0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0  };
    static const __m128i xmm_increment_j = { 0x0000000400000004ull, 0x0000000400000004ull };

    static const __m256d ymm_merge_mask_table[4] = {
        { +0.0, +0.0, +0.0, +0.0 },
        { -0.0, -0.0, -0.0, +0.0 },
        { -0.0, -0.0, +0.0, +0.0 },
        { -0.0, +0.0, +0.0, +0.0 }
    };
#endif

#ifdef __MIC__
    ERD_OFFLOAD static const __m512d zmm_c0 = {  0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8  };
    ERD_OFFLOAD static const __m512d zmm_c1 = {  0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0  };
    ERD_OFFLOAD static const __m512d zmm_c2 = { -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1  };
    ERD_OFFLOAD static const __m512d zmm_c3 = {  0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3  };
    ERD_OFFLOAD static const __m512d zmm_c4 = { -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6  };
    ERD_OFFLOAD static const __m512d zmm_c5 = {  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10 };
    ERD_OFFLOAD static const __m512d zmm_x0 = {  0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5  };
    ERD_OFFLOAD static const __m512d zmm_one = { 0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0  };

    ERD_OFFLOAD static const __m512i zmm_increment_j = { 8u, 8u, 8u, 8u, 8u, 8u, 8u, 8u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u };
    ERD_OFFLOAD static const __m512i zmm_increment_i = { 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u };
    ERD_OFFLOAD static const __m512i zmm_init_j = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif

ERD_OFFLOAD __attribute__((noinline)) uint32_t set_pairs(
    uint32_t npgtoa, uint32_t npgtob, double rnabsq,
    const double alphaa[restrict static npgtoa], const double alphab[restrict static npgtoa],
    uint32_t prima[restrict static npgtoa*npgtob], uint32_t primb[restrict static npgtoa*npgtob],
    double rho[restrict static npgtoa*npgtob],
    double qmin, double smaxcd, double rminsq)
{
    const double csmaxcd = smaxcd * (0x1.C5BF891B4EF6Bp-1 / TOL);
    uint32_t nij = 0;
    if (npgtoa > npgtob) {
        const double *restrict alphaa_copy = alphaa;
        alphaa = alphab;
        alphab = alphaa_copy;
        uint32_t *restrict prima_copy = prima;
        prima = primb;
        primb = prima_copy;
        uint32_t npgtoa_copy = npgtoa;
        npgtoa = npgtob;
        npgtob = npgtoa_copy;
    }
    const double minus_rnabsq = -rnabsq;
    const double rminsq_qmin = rminsq * qmin;
    for (uint32_t i = 0; i < npgtoa; i += 1) {
        const double a = alphaa[i];
        for (uint32_t j = 0; j < npgtob; j += 1) {
            const double b = alphab[j];
            const double p = a + b;
            const double ab = a * b;
            const double pqp = p + qmin;
            const double pqpinv = 1.0 / pqp;
            const double rhoab = __builtin_exp(ab * minus_rnabsq / p);
            const double t = rminsq_qmin * (p / pqp);

            /* pi/4 * f0(x) == x */
            const double x = (t == 0.0) ? x0 : t;
            /* approximates square(erf(sqrt(x))) on [0, 5] */
            const double f0 = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))));
            const double f = f0 < 1.0 ? f0 : 1.0;
            if (1){//ab*square(square(rhoab*csmaxcd) * (ab*f)) >= square((x*pqp)*square(p))) {
                rho[nij] = rhoab;
                prima[nij] = i;
                primb[nij] = j;
                nij += 1;
            }
        }
    }
    return nij;
}

ERD_OFFLOAD void erd__set_ij_kl_pairs(
    uint32_t npgtoa, uint32_t npgtob, uint32_t npgtoc, uint32_t npgtod,
    double minalphaa, double minalphab, double minalphac, double minalphad,
    double xa, double ya, double za,
    double xb, double yb, double zb,
    double xc, double yc, double zc,
    double xd, double yd, double zd,
    double rnabsq, double rncdsq, double prefact,
    const double alphaa[restrict static npgtoa],
    const double alphab[restrict static npgtob],
    const double alphac[restrict static npgtoc],
    const double alphad[restrict static npgtod],
    uint32_t nij_ptr[restrict static 1], uint32_t nkl_ptr[restrict static 1],
    uint32_t prima[restrict static npgtoa*npgtob], uint32_t primb[restrict static npgtoa*npgtob], uint32_t primc[restrict static npgtoc*npgtod], uint32_t primd[restrict static npgtoc*npgtod],
    double rhoab[restrict static npgtoa*npgtob],
    double rhocd[restrict static npgtoc*npgtod])
{
    // compute min
    const double rminsq = erd__dsqmin_line_segments(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd);

    const double pmin = minalphaa + minalphab;
    const double abmin = minalphaa * minalphab;
    const double pinv = 1.0 / pmin;
    const double smaxab = prefact * pow3o4(abmin) * __builtin_exp(-abmin * rnabsq * pinv) * pinv;

    const double qmin = minalphac + minalphad;
    const double cdmin = minalphac * minalphad;
    const double qinv = 1.0 / qmin;
    const double smaxcd = prefact * pow3o4(cdmin) * __builtin_exp(-cdmin * rncdsq * qinv) * qinv;

    /* ...perform K2 primitive screening on A,B part. */
    uint32_t nij = set_pairs(npgtoa, npgtob, rnabsq, alphaa, alphab, prima, primb, rhoab, qmin, smaxcd, rminsq);
    if (nij == 0) {
        *nij_ptr = 0;
        *nkl_ptr = 0;
        return;
    }

    uint32_t nkl = set_pairs(npgtoc, npgtod, rncdsq, alphac, alphad, primc, primd, rhocd, pmin, smaxab, rminsq);
    if (nkl == 0) {
        *nij_ptr = 0;
        *nkl_ptr = 0;
        return;
    }

    *nij_ptr = nij;
    *nkl_ptr = nkl;
}
