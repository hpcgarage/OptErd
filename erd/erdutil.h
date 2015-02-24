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
#include <stdint.h>
#include <stdbool.h>

#ifdef __x86_64__
    #include <x86intrin.h>
#elif defined(__MIC__)
    #include <immintrin.h>
#endif

#ifdef __INTEL_OFFLOAD
    #define ERD_OFFLOAD __attribute__((target(mic)))
#else
    #define ERD_OFFLOAD
#endif

#ifdef __ERD_PROFILE__
    #define ERD_PROFILE_START(function) \
        const uint64_t start_ticks_##function = __rdtsc();
    #define ERD_PROFILE_END(function) \
        const uint64_t end_ticks_##function = __rdtsc(); \
        erd_ticks[tid][function##_ticks] += end_ticks_##function - start_ticks_##function;
#else
    #define ERD_PROFILE_START(function)
    #define ERD_PROFILE_END(function)
#endif

#if defined(__MIC__) || defined(__AVX512F__)
    #define ERD_SIMD_SIZE 64
    #define ERD_SIMD_WIDTH 8
#elif defined(__AVX__)
    #define ERD_SIMD_SIZE 32
    #define ERD_SIMD_WIDTH 4
#elif defined(__SSE2__)
    #define ERD_SIMD_SIZE 16
    #define ERD_SIMD_WIDTH 2
#else
    #define ERD_SIMD_SIZE 8
    #define ERD_SIMD_WIDTH 1
#endif

#define ERD_ALIGN(alignment) __attribute__((aligned(alignment)))
#define ERD_SIMD_ALIGN ERD_ALIGN(ERD_SIMD_SIZE)

#ifdef __INTEL_COMPILER
    #define ERD_ASSUME_ALIGNED(pointer, alignment) __assume_aligned(pointer, alignment);
#else
    #define ERD_ASSUME_ALIGNED(pointer, alignment) pointer = __builtin_assume_aligned(pointer, alignment);
#endif

#if defined(__MIC__) || defined(__AVX512F__)
    #define ERD_SIMD_ZERO_TAIL_64f(array, simd_length) _mm512_store_pd(&array[simd_length - ERD_SIMD_WIDTH], _mm512_setzero_pd())
#elif defined(__AVX__)
    #define ERD_SIMD_ZERO_TAIL_64f(array, simd_length) _mm256_store_pd(&array[simd_length - ERD_SIMD_WIDTH], _mm256_setzero_pd())
#elif defined(__SSE2__)
    #define ERD_SIMD_ZERO_TAIL_64f(array, simd_length) _mm_store_pd(&array[simd_length - ERD_SIMD_WIDTH], _mm_setzero_pd())
#else
    #define ERD_SIMD_ZERO_TAIL_64f(array, simd_length)
#endif

#define ERD_SWAP(x, y) \
    ({ __typeof__(x) __temp = x; \
    x = y; \
    y = __temp; })

ERD_OFFLOAD static inline uint32_t max32u(uint32_t a, uint32_t b) {
    return a > b ? a : b;
}

ERD_OFFLOAD static inline uint32_t max3x32u(uint32_t a, uint32_t b, uint32_t c) {
    return max32u(max32u(a, b), c);
}

ERD_OFFLOAD static inline uint32_t max4x32u(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    return max32u(max32u(a, b), max32u(c, d));
}

ERD_OFFLOAD static inline uint32_t min32u(uint32_t a, uint32_t b) {
    return a <= b ? a : b;
}

ERD_OFFLOAD static inline uint32_t min3x32u(uint32_t a, uint32_t b, uint32_t c) {
    return min32u(min32u(a, b), c);
}

ERD_OFFLOAD static inline uint32_t min4x32u(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    return min32u(min32u(a, b), min32u(c, d));
}
