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

#pragma once

#include <stddef.h>

#define NGRID    920
#define MGRID     10

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

static const double tmax = 46.0;
static const double tvstep = 20.0;
static const double tstep = 0.05;

extern double boys_table[NGRID + 1][MGRID + 1];

struct Boys01 {
    double f0;
    double f1;
};

struct Boys012 {
    double f0;
    double f1;
    double f2;
};

struct Boys0123 {
    double f0;
    double f1;
    double f2;
    double f3;
};

struct Boys01234 {
    double f0;
    double f1;
    double f2;
    double f3;
    double f4;
};

static inline double boys0(double t, double scale) {
    double f0;
    if (t <= tmax) {
        const int tgrid = __builtin_lround(t * tvstep);
        const double delta1 = tgrid * tstep - t;
        const double delta2 = delta1 * 0x1.0000000000000p-1;
        const double delta3 = delta1 * 0x1.5555555555555p-2;
        const double delta4 = delta1 * 0x1.0000000000000p-2;
        const double delta5 = delta1 * 0x1.999999999999Ap-3;
        const double delta6 = delta1 * 0x1.5555555555555p-3;
        f0 = (((((boys_table[tgrid][6] * delta6 +
            boys_table[tgrid][5]) * delta5 +
                boys_table[tgrid][4]) * delta4 +
                    boys_table[tgrid][3]) * delta3 +
                        boys_table[tgrid][2]) * delta2 +
                            boys_table[tgrid][1]) * delta1 +
                                boys_table[tgrid][0];
    } else {
        /* sqrt(pi) / 2 */
        const double factor = 0x1.C5BF891B4EF6Bp-1;
        const double tinv = 1.0 / t;
        f0 = factor * __builtin_sqrt(tinv);
    }
    return scale * f0;
}

static inline struct Boys01 boys01(double t, double scale) {
    struct Boys01 boys;
    if (t <= tmax) {
        const int tgrid = __builtin_lround(t * tvstep);
        const double delta1 = tgrid * tstep - t;
        const double delta2 = delta1 * 0x1.0000000000000p-1;
        const double delta3 = delta1 * 0x1.5555555555555p-2;
        const double delta4 = delta1 * 0x1.0000000000000p-2;
        const double delta5 = delta1 * 0x1.999999999999Ap-3;
        const double delta6 = delta1 * 0x1.5555555555555p-3;
        const double f0 = (((((boys_table[tgrid][6] * delta6 +
            boys_table[tgrid][5]) * delta5 +
                boys_table[tgrid][4]) * delta4 +
                    boys_table[tgrid][3]) * delta3 +
                        boys_table[tgrid][2]) * delta2 +
                            boys_table[tgrid][1]) * delta1 +
                                boys_table[tgrid][0];
        const double f1 = (((((boys_table[tgrid][7] * delta6 +
            boys_table[tgrid][6]) * delta5 +
                boys_table[tgrid][5]) * delta4 +
                    boys_table[tgrid][4]) * delta3 +
                        boys_table[tgrid][3]) * delta2 +
                            boys_table[tgrid][2]) * delta1 +
                                boys_table[tgrid][1];
        boys.f0 = scale * f0;
        boys.f1 = scale * f1;
    } else {
        /* sqrt(pi) / 2 */
        const double factor = 0x1.C5BF891B4EF6Bp-1;
        const double tinv = 1.0 / t;
        boys.f0 = (scale * factor) * __builtin_sqrt(tinv);
        boys.f1 = (tinv * 0.5) * boys.f0;
    }
    return boys;
}

static inline struct Boys012 boys012(double t, double scale) {
    struct Boys012 boys;
    if (t <= tmax) {
        const int tgrid = __builtin_lround(t * tvstep);
        const double delta1 = tgrid * tstep - t;
        const double delta2 = delta1 * 0x1.0000000000000p-1;
        const double delta3 = delta1 * 0x1.5555555555555p-2;
        const double delta4 = delta1 * 0x1.0000000000000p-2;
        const double delta5 = delta1 * 0x1.999999999999Ap-3;
        const double delta6 = delta1 * 0x1.5555555555555p-3;
        const double f0 = (((((boys_table[tgrid][6] * delta6 +
            boys_table[tgrid][5]) * delta5 +
                boys_table[tgrid][4]) * delta4 +
                    boys_table[tgrid][3]) * delta3 +
                        boys_table[tgrid][2]) * delta2 +
                            boys_table[tgrid][1]) * delta1 +
                                boys_table[tgrid][0];
        const double f1 = (((((boys_table[tgrid][7] * delta6 +
            boys_table[tgrid][6]) * delta5 +
                boys_table[tgrid][5]) * delta4 +
                    boys_table[tgrid][4]) * delta3 +
                        boys_table[tgrid][3]) * delta2 +
                            boys_table[tgrid][2]) * delta1 +
                                boys_table[tgrid][1];
        const double f2 = (((((boys_table[tgrid][8] * delta6 +
            boys_table[tgrid][7]) * delta5 +
                boys_table[tgrid][6]) * delta4 +
                    boys_table[tgrid][5]) * delta3 +
                        boys_table[tgrid][4]) * delta2 +
                            boys_table[tgrid][3]) * delta1 +
                                boys_table[tgrid][2];
        boys.f0 = scale * f0;
        boys.f1 = scale * f1;
        boys.f2 = scale * f2;
    } else {
        /* sqrt(pi) / 2 */
        const double factor = 0x1.C5BF891B4EF6Bp-1;
        const double tinv = 1.0 / t;
        boys.f0 = (scale * factor) * __builtin_sqrt(tinv);
        boys.f1 = (tinv * 0.5) * boys.f0;
        boys.f2 = (tinv * 1.5) * boys.f1;
    }
    return boys;
}

static inline struct Boys0123 boys0123(double t, double scale) {
    struct Boys0123 boys;
    if (t <= tmax) {
        const int tgrid = __builtin_lround(t * tvstep);
        const double delta1 = tgrid * tstep - t;
        const double delta2 = delta1 * 0x1.0000000000000p-1;
        const double delta3 = delta1 * 0x1.5555555555555p-2;
        const double delta4 = delta1 * 0x1.0000000000000p-2;
        const double delta5 = delta1 * 0x1.999999999999Ap-3;
        const double delta6 = delta1 * 0x1.5555555555555p-3;
        const double f0 = (((((boys_table[tgrid][6] * delta6 +
            boys_table[tgrid][5]) * delta5 +
                boys_table[tgrid][4]) * delta4 +
                    boys_table[tgrid][3]) * delta3 +
                        boys_table[tgrid][2]) * delta2 +
                            boys_table[tgrid][1]) * delta1 +
                                boys_table[tgrid][0];
        const double f1 = (((((boys_table[tgrid][7] * delta6 +
            boys_table[tgrid][6]) * delta5 +
                boys_table[tgrid][5]) * delta4 +
                    boys_table[tgrid][4]) * delta3 +
                        boys_table[tgrid][3]) * delta2 +
                            boys_table[tgrid][2]) * delta1 +
                                boys_table[tgrid][1];
        const double f2 = (((((boys_table[tgrid][8] * delta6 +
            boys_table[tgrid][7]) * delta5 +
                boys_table[tgrid][6]) * delta4 +
                    boys_table[tgrid][5]) * delta3 +
                        boys_table[tgrid][4]) * delta2 +
                            boys_table[tgrid][3]) * delta1 +
                                boys_table[tgrid][2];
        const double f3 = (((((boys_table[tgrid][9] * delta6 +
            boys_table[tgrid][8]) * delta5 +
                boys_table[tgrid][7]) * delta4 +
                    boys_table[tgrid][6]) * delta3 +
                        boys_table[tgrid][5]) * delta2 +
                            boys_table[tgrid][4]) * delta1 +
                                boys_table[tgrid][3];
        boys.f0 = scale * f0;
        boys.f1 = scale * f1;
        boys.f2 = scale * f2;
        boys.f3 = scale * f3;
    } else {
        /* sqrt(pi) / 2 */
        const double factor = 0x1.C5BF891B4EF6Bp-1;
        const double tinv = 1.0 / t;
        boys.f0 = (scale * factor) * __builtin_sqrt(tinv);
        boys.f1 = (tinv * 0.5) * boys.f0;
        boys.f2 = (tinv * 1.5) * boys.f1;
        boys.f3 = (tinv * 2.5) * boys.f2;
    }
    return boys;
}

static inline struct Boys01234 boys01234(double t, double scale) {
    struct Boys01234 boys;
    if (t <= tmax) {
        const int tgrid = __builtin_lround(t * tvstep);
        const double delta1 = tgrid * tstep - t;
        const double delta2 = delta1 * 0x1.0000000000000p-1;
        const double delta3 = delta1 * 0x1.5555555555555p-2;
        const double delta4 = delta1 * 0x1.0000000000000p-2;
        const double delta5 = delta1 * 0x1.999999999999Ap-3;
        const double delta6 = delta1 * 0x1.5555555555555p-3;
        const double f0 = (((((boys_table[tgrid][6] * delta6 +
            boys_table[tgrid][5]) * delta5 +
                boys_table[tgrid][4]) * delta4 +
                    boys_table[tgrid][3]) * delta3 +
                        boys_table[tgrid][2]) * delta2 +
                            boys_table[tgrid][1]) * delta1 +
                                boys_table[tgrid][0];
        const double f1 = (((((boys_table[tgrid][7] * delta6 +
            boys_table[tgrid][6]) * delta5 +
                boys_table[tgrid][5]) * delta4 +
                    boys_table[tgrid][4]) * delta3 +
                        boys_table[tgrid][3]) * delta2 +
                            boys_table[tgrid][2]) * delta1 +
                                boys_table[tgrid][1];
        const double f2 = (((((boys_table[tgrid][8] * delta6 +
            boys_table[tgrid][7]) * delta5 +
                boys_table[tgrid][6]) * delta4 +
                    boys_table[tgrid][5]) * delta3 +
                        boys_table[tgrid][4]) * delta2 +
                            boys_table[tgrid][3]) * delta1 +
                                boys_table[tgrid][2];
        const double f3 = (((((boys_table[tgrid][9] * delta6 +
            boys_table[tgrid][8]) * delta5 +
                boys_table[tgrid][7]) * delta4 +
                    boys_table[tgrid][6]) * delta3 +
                        boys_table[tgrid][5]) * delta2 +
                            boys_table[tgrid][4]) * delta1 +
                                boys_table[tgrid][3];
        const double f4 = (((((boys_table[tgrid][10] * delta6 +
            boys_table[tgrid][9]) * delta5 +
                boys_table[tgrid][8]) * delta4 +
                    boys_table[tgrid][7]) * delta3 +
                        boys_table[tgrid][6]) * delta2 +
                            boys_table[tgrid][5]) * delta1 +
                                boys_table[tgrid][4];
        boys.f0 = scale * f0;
        boys.f1 = scale * f1;
        boys.f2 = scale * f2;
        boys.f3 = scale * f3;
        boys.f4 = scale * f4;
    } else {
        /* sqrt(pi) / 2 */
        const double factor = 0x1.C5BF891B4EF6Bp-1;
        const double tinv = 1.0 / t;
        boys.f0 = (scale * factor) * __builtin_sqrt(tinv);
        boys.f1 = (tinv * 0.5) * boys.f0;
        boys.f2 = (tinv * 1.5) * boys.f1;
        boys.f3 = (tinv * 2.5) * boys.f2;
        boys.f4 = (tinv * 3.5) * boys.f3;
    }
    return boys;
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
