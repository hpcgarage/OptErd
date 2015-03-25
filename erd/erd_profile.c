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
#include <unistd.h>

#include "erd_profile.h"
#include "erd.h"
#include "erdutil.h"


ERD_CACHE_ALIGN uint64_t erd_ticks[MAXTHREADS][erd__num_ticks + ERD_CACHELINE/sizeof(uint64_t)];
__thread int erd_tid = 0;

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


void erd_reset_profile(int nthreads)
{
    printf("Reset profiler for %d threads\n", nthreads);
    for (int tid = 0; tid < nthreads; tid++) {
        memset(erd_ticks[tid], 0, sizeof(uint64_t) * erd__num_ticks);
    }
}


void erd_print_profile(int nthreads, int mode)
{
    double total_secs = 0.0;
    const uint64_t start_clock = ReadTSC();
    sleep(1);
    const uint64_t end_clock = ReadTSC();
    const uint64_t freq = end_clock - start_clock;

    printf("\n");    
    printf("ticks = %.3e\n", (double)freq/1e9);
    
    for (int k = 0; k < erd__num_ticks; k++) {
        total_secs = 0.0;
        printf("%28s", ticks_name[k]);
        for (int i = 0; i < nthreads; i++) {
            if (mode == 0) {
                printf("\t%.3f", (double)erd_ticks[i][k]/freq);
            }
            total_secs += (double)erd_ticks[i][k]/freq;
        }        
        printf(",\t%.3f", total_secs/nthreads);
        printf("\n");
    }
}