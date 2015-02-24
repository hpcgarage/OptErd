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

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>

#include "erd_profile.h"


__declspec(align(256)) uint64_t erd_ticks[MAXTHREADS][erd__num_ticks + 8];

static char ticks_name[erd__num_ticks][128] = 
{
    "erd__prepare_ctr",
    "erd__set_ij_kl_pairs",
    "erd__rys_roots_weights",
    "erd__2d_coefficients",
    "erd__2d_pq_integrals",
    "erd__int2d_to_e0f0",
    "@erd__e0f0_pcgto_block",
    "erd__ctr_4index_block",
    "erd__xyz_to_ry_abcd",
    "erd__hrr_matrix",
    "erd__hrr_transform",
    "erd__move_ry",
    "erd__spherical_transform",
    "@erd__csgto",

    "erd__1111_prepare_ctr",
    "erd__1111_set_ij_kl_pairs",    
    "erd__ssss_pcgto_block",
    "erd__sssp_pcgto_block",
    "erd__sspp_pcgto_block",
    "erd__sppp_pcgto_block",
    "erd__pppp_pcgto_block",
    "erd__1111_ctr_4index_block",
    "@erd__1111_csgto"
};


void erd_reset_profile(void) {
    #ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
    assert (nthreads <= MAXTHREADS);
    #else
    const int nthreads = 1;
    #endif
    printf("Reset profiler for %d threads\n", nthreads);
    #pragma omp parallel for
    for (int i = 0; i < nthreads; i++) {
        memset (erd_ticks[i], 0, sizeof(uint64_t) * erd__num_ticks);
    }
}

void erd_print_profile(int mode) {
    double total_secs = 0.0;

    #ifdef _OPENMP
    const int nthreads = omp_get_max_threads ();
    #else
    const int nthreads = 1;
    #endif
    const uint64_t start_clock = __rdtsc();
    sleep(1);
    const uint64_t end_clock = __rdtsc();
    const uint64_t freq = end_clock - start_clock;

    printf ("\n");    
    printf("freq = %.3lf GHz\n", (double)freq/1e9); 
    for (int k = 0; k < erd__num_ticks; k++) {
        total_secs = 0.0;
        printf("%28s", ticks_name[k]);
        for (int i = 0; i < nthreads; i++) {
            if (mode == 0) {
                printf("\t%.3lf", (double)erd_ticks[i][k]/freq);
            }
            total_secs += (double)erd_ticks[i][k]/freq;
        }        
        printf(",\t%.3lf", total_secs/nthreads);
        printf("\n");
    }
}


int start_timer_ (uint64_t *stime) {
    *stime  = __rdtsc();

    return 0;
}


int end_timer_ (int *idx, uint64_t *stime)
{
    #ifdef _OPENMP
    const int tid = omp_get_thread_num();   
    #else
    const int tid = 0;
    #endif
    const uint64_t etime = __rdtsc();
    erd_ticks[tid][*idx] += etime - *stime;

    return 0;
}
