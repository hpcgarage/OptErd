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

#ifndef ERD_PROFILE_H_
#define ERD_PROFILE_H_


#include <stdint.h>


#define MAXTHREADS     512


typedef enum
{
    // general case
    erd__prepare_ctr_ticks              = 0,
    erd__set_ij_kl_pairs_ticks          = 1,
    // VRR
    erd__rys_roots_weights_ticks        = 2,
    erd__2d_coefficients_ticks          = 3,
    erd__2d_pq_integrals_ticks          = 4,
    erd__int2d_to_e0f0_ticks            = 5,
    erd__e0f0_pcgto_block_ticks         = 6,
    erd__ctr_4index_block_ticks         = 7,
    
    erd__xyz_to_ry_abcd_ticks           = 8,
    erd__hrr_matrix_ticks               = 9,
    erd__hrr_transform_ticks            = 10,
    erd__move_ry_ticks                  = 11,
    erd__spherical_transform_ticks      = 12,
    erd__csgto_ticks                    = 13,

    // 1111 case
    erd__1111_prepare_ctr_ticks         = 14,
    erd__1111_set_ij_kl_pairs_ticks     = 15,
    erd__ssss_pcgto_block_ticks         = 16,
    erd__sssp_pcgto_block_ticks         = 17,
    erd__sspp_pcgto_block_ticks         = 18,
    erd__sppp_pcgto_block_ticks         = 19,
    erd__pppp_pcgto_block_ticks         = 20,
    erd__1111_ctr_4index_block_ticks    = 21,
    erd__1111_csgto_ticks               = 22,
    erd__num_ticks
} ErdTicks_t;


extern uint64_t erd_ticks[MAXTHREADS][erd__num_ticks + 8];

void erd_reset_profile(void);

void erd_print_profile(int mode);


inline uint64_t ReadTSC(void)
{
#if defined(__i386__)
    uint64_t x;
    __asm__ __volatile__ (".byte 0x0f, 0x31" : "=A" (x));
    return x;
#elif defined(__x86_64__)
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (uint64_t)lo)|( ((uint64_t)hi)<<32 );
#elif defined(__powerpc__)
    uint64_t result=0;
    uint64_t upper, lower,tmp;
    __asm__ __volatile__ (
        "0:                  \n"
        "\tmftbu   %0           \n"
        "\tmftb    %1           \n"
        "\tmftbu   %2           \n"
        "\tcmpw    %2,%0        \n"
        "\tbne     0b         \n"
        : "=r"(upper),"=r"(lower),"=r"(tmp)
    );
    result = upper;
    result = result<<32;
    result = result|lower;
    return result;
#else
    return 0;
#endif // defined(__i386__)
}


#if defined(ERD_PROFILE_) && defined(__linux__)
    #define ERD_PROFILE_START(function) \
    #ifdef _OPENMP \
        {const int tid = omp_get_thread_num(); \
    #else \
        const int tid = 0; \
    #endif \
        const uint64_t start_ticks_##function = ReadTSC();
    #define ERD_PROFILE_END(function) \
        const uint64_t end_ticks_##function = ReadTSC(); \
        erd_ticks[tid][function##_ticks] += end_ticks_##function - start_ticks_##function;}
#else
    #define ERD_PROFILE_START(function)
    #define ERD_PROFILE_END(function)
#endif


#endif // ERD_PROFILE_H_
