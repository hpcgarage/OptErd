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

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

/* (2i+1) / ((4i-1) * (4i+3)) for i = 1...100 */
extern const double r2[100];

/* ((4i-3) * (4i+1) * square(4i-1)) / (2i * (2i+1) * square(2i-1)) for i = 1...100 */
extern const double sinv[100];

extern const double csmall[16];

/* (4i * (2i+1) - 1) / ((4i+3) * (4i-1)) for i = 0...99 */
extern const double ajac[100];

/* (4*square(i) * (4*square(i) - 4*i + 1)) / ((4i-3) * (4i+1) * square(4i-1)) for i = 1...99 */
extern const double bjac[99];

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
