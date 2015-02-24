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

#ifndef __OED_INTEGRAL_H__
#define __OED_INTEGRAL_H__


#define OED_SPHERIC 1
#define OED_SCREEN 0


extern void oed__gener_kin_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2,
                                   int *npgto1, int *npgto2,
                                   int *shell1, int *shell2,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   double *alpha, double *cc,
                                   int *ccbeg, int *ccend, int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst, double *zcore);


extern void oed__gener_nai_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2,
                                   int *npgto1, int *npgto2,
                                   int *shell1, int *shell2,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   int *natoms,
                                   double *xn, double *yn, double *zn,
                                   double *charge, double *alpha, double *cc,
                                   int *ccbeg, int *ccend,
                                   int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst, double *zcore);


extern void oed__gener_ovl_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2,
                                   int *npgto1, int *npgto2,
                                   int *shell1, int *shell2,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   double *alpha, double *cc, int *ccbeg, int *ccend,
                                   int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst, double *zcore);


extern void oed__memory_ovl_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2,
                                    int *npgto1, int *npgto2,
                                    int *shell1, int *shell2,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt,
                                    int *zmin, int *zopt);


extern void oed__memory_kin_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2,
                                    int *npgto1, int *npgto2,
                                    int *shell1, int *shell2,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt, int *zmin, int *zopt);

extern void oed__memory_nai_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2,
                                    int *npgto1, int *npgto2,
                                    int *shell1, int *shell2,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    int *natoms, double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt, int *zmin, int *zopt);


#endif /* __OED_INTEGRAL_H__ */
