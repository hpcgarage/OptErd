/*
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

#define ERD_SCREEN true
#define ERD_SPHERIC 1
#define ERD_CARTESIAN 0


#define MAX(a,b)    ((a) < (b) ? (b) : (a))


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

extern uint32_t erd__1111_csgto(
    uint32_t A, uint32_t B, uint32_t C, uint32_t D,
    const uint32_t npgto[restrict static 1], const uint32_t shell[restrict static 1], const double xyz0[restrict static 1],
    const double *restrict alpha[restrict static 1], const double minalpha[restrict static 1], const double *restrict cc[restrict static 1], const double *restrict norm[restrict static 1],
    uint32_t buffer_capacity, uint32_t integral_counts[restrict static 1], double output_buffer[restrict static 81]);

extern uint32_t erd__csgto(
    uint32_t A, uint32_t B, uint32_t C, uint32_t D,
    const uint32_t npgto[restrict static 1], const uint32_t shell[restrict static 1], const double xyz0[restrict static 1],
    const double *restrict alpha[restrict static 1], const double minalpha[restrict static 1], const double *restrict cc[restrict static 1], const double *restrict norm[restrict static 1],
    int **vrrtab,
    bool spheric,
    uint32_t buffer_capacity, uint32_t integral_counts[restrict static 1], double output_buffer[restrict static 1]);

extern size_t erd__memory_csgto(uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    bool spheric);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
