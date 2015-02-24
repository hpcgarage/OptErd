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

#ifndef __CINT_BASISSET_H__
#define __CINT_BASISSET_H__


#include <cint_def.h>

void _maxMomentum (BasisSet_t basis, int *max_momentum);

void _maxPrimid (BasisSet_t basis, int *max_primid);

void _maxnumExp (BasisSet_t basis, int *max_nexp);

CIntStatus_t import_basis (char *file, BasisSet_t basis);

CIntStatus_t import_molecule (char *file, BasisSet_t basis);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

CIntStatus_t parse_molecule (BasisSet_t basis);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#endif /* __CINT_BASISSET_H__ */
