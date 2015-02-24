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

#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include "boys.h"
#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SSSP_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between s-shell and p-shell primitive */
/*                spherical gaussian type orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals, that is its dimension is determined by */
/*                the total number of primitive functions (here = 3) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                        sssp , ssps , spss , psss */
/*                The cartesian primitive integrals are evaluated using */
/*                the auxillary functions technique, described by */
/*                V.R. Saunders, "An Introduction to Molecular Integral */
/*                Evaluation" in "Computational Techniques in Quantum */
/*                Chemistry and Molecular Physics" edited by */
/*                GHF Diercksen, BT Sutcliff and A Veillard, D. Reidel */
/*                Publ. Comp. Dordrecht (1975), p. 347. */
/*                The cartesian primitive integrals are each evaluated */
/*                explicitely using common multipliers if possible and */
/*                avoiding array addresses. */
/*                  Input: */
/*                    NBATCH       =  size of the primitive cartesian */
/*                                    sssp/ssps/spss/psss integral batch */
/*                    ATOMIC       =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                    ATOM12(34)   =  indicates, if centers 1 and 2 */
/*                                    (3 and 4) coincide */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the contracted shell pairs 1,2 */
/*                                    (3,4) */
/*                    NIJ          =  total # of ij primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair 1,2 */
/*                    NIJBEG(END)  =  first(last) ij primitive index */
/*                                    defining the ij block */
/*                    NKL          =  total # of kl primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair 3,4 */
/*                    NKLBEG(END)  =  first(last) kl primitive index */
/*                                    defining the kl block */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = 1,2,3,4 */
/*                    SHELLx       =  the shell type for contraction */
/*                                    shells x = 1,3,P=1+2 */
/*                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers */
/*                                    x = 1,2,3,4 */
/*                    Xxx,Yxx,Zxx  =  the x,y,z-coordinate differences */
/*                                    between centers xx = 12 and 34 */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = 1,2,3,4 */
/*                    FTABLE       =  Fm (T) table for interpolation */
/*                                    in low T region */
/*                    MGRID        =  maximum m in Fm (T) table */
/*                    NGRID        =  # of T's for which Fm (T) table */
/*                                    was set up */
/*                    TMAX         =  maximum T in Fm (T) table */
/*                    TSTEP        =  difference between two consecutive */
/*                                    T's in Fm (T) table */
/*                    TVSTEP       =  Inverse of TSTEP */
/*                    PRIMx        =  i,j,k,l labels of primitives for */
/*                                    the respective contraction shells */
/*                                    x = 1,2,3,4 */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = 1,2,3,4 */
/*                    RHO12(34)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells 1 and 2 */
/*                                    (3 and 4) */
/*                    P            =  will hold current MIJ exponent */
/*                                    sums for contraction shells 1 */
/*                                    and 2 */
/*                    Px           =  will hold current MIJ coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers P=1+2 */
/*                    SCALEP       =  will hold current MIJ values of */
/*                                    scaling factors related to point P */
/*                    Q            =  will hold current MKL exponent */
/*                                    sums for contraction shells 3 */
/*                                    and 4 */
/*                    Qx           =  will hold current MKL coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers Q=3+4 */
/*                    SCALEQ       =  will hold current MKL values of */
/*                                    scaling factors related to point Q */
/*                  Output: */
/*                    BATCH        =  current batch of primitive */
/*                                    cartesian sssp/ssps/spss/psss */
/*                                    integrals */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__sssp_pcgto_block(uint32_t nij, uint32_t nkl,
    uint32_t shell1, uint32_t shell3, uint32_t shellp,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double *restrict alpha1, const double *restrict alpha2, const double *restrict alpha3, const double *restrict alpha4,
    const double *restrict cc1, const double *restrict cc2, const double *restrict cc3, const double *restrict cc4,
    const uint32_t *restrict prim1, const uint32_t *restrict prim2, const uint32_t *restrict prim3, const uint32_t *restrict prim4,
    const double *restrict norm1, const double *restrict norm2, const double *restrict norm3, const double *restrict norm4,
    const double *restrict rho12, const double *restrict rho34,
    double *restrict cbatch)
{
    const double x12 = x1 - x2;
    const double y12 = y1 - y2;
    const double z12 = z1 - z2;
    const double x34 = x3 - x4;
    const double y34 = y3 - y4;
    const double z34 = z3 - z4;

    ERD_SIMD_ALIGN double p[PAD_LEN(nij)], px[PAD_LEN(nij)], py[PAD_LEN(nij)], pz[PAD_LEN(nij)], scalep[PAD_LEN(nij)];
    for (uint32_t ij = 0; ij < nij; ij += 1) {
        const uint32_t i = prim1[ij];
        const uint32_t j = prim2[ij];
        const double exp1 = alpha1[i];
        const double exp2 = alpha2[j];
        double pval = exp1 + exp2;
        p[ij] = pval;
        pval = exp1 / pval;
        px[ij] = pval * x12 + x2;
        py[ij] = pval * y12 + y2;
        pz[ij] = pval * z12 + z2;
        scalep[ij] = cc1[i] * cc2[j] * norm1[i] * norm2[j] * rho12[ij];
    }

    ERD_SIMD_ALIGN double q[PAD_LEN(nkl)], qx[PAD_LEN(nkl)], qy[PAD_LEN(nkl)], qz[PAD_LEN(nkl)], scaleq[PAD_LEN(nkl)];
    for (uint32_t kl = 0; kl < nkl; kl += 1) {
        const uint32_t k = prim3[kl];
        const uint32_t l = prim4[kl];
        const double exp3 = alpha3[k];
        const double exp4 = alpha4[l];
        double qval = exp3 + exp4;
        q[kl] = qval;
        qval = exp3 / qval;
        qx[kl] = qval * x34 + x4;
        qy[kl] = qval * y34 + y4;
        qz[kl] = qval * z34 + z4;
        scaleq[kl] = cc3[k] * cc4[l] * norm3[k] * norm4[l] * rho34[kl];
    }

    if (shellp == 0) {
        // 0     5   |  (AB|CD)  4-center   sssp and ssps
        double qxsub, qysub, qzsub;
        if (shell3 == 1) {
            qxsub = x3;
            qysub = y3;
            qzsub = z3;
        } else {
            qxsub = x4;
            qysub = y4;
            qzsub = z4;
        }
        for (uint32_t ij = 0; ij < nij; ij += 1) {
            const double pval = p[ij];
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double pscale = scalep[ij];
            for (uint32_t kl = 0; kl < nkl; kl += 1) {
                const double qval = q[kl];
                const double qxval = qx[kl];
                const double qyval = qy[kl];
                const double qzval = qz[kl];
                const double pqmult = pval * qval;
                const double pqplus = pval + qval;
                const double pqpinv = 1. / pqplus;
                const double pqx = pxval - qxval;
                const double pqy = pyval - qyval;
                const double pqz = pzval - qzval;
                const double t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
                const double scale = pscale * scaleq[kl] / (pqmult * __builtin_sqrt(pqplus));
                struct Boys01 boys = boys01(t, scale);
                boys.f1 = boys.f1 * pval * pqpinv;
                cbatch[0] += (qxval - qxsub) * boys.f0 + pqx * boys.f1;
                cbatch[1] += (qyval - qysub) * boys.f0 + pqy * boys.f1;
                cbatch[2] += (qzval - qzsub) * boys.f0 + pqz * boys.f1;
            }
        }
    } else {
        // 1     5   |  (AB|CD)  4-center   spss and psss
        double pxsub, pysub, pzsub;
        if (shell1 == 1) {
            pxsub = x1;
            pysub = y1;
            pzsub = z1;
        } else {
            pxsub = x2;
            pysub = y2;
            pzsub = z2;
        }
        for (uint32_t ij = 0; ij < nij; ij += 1) {
            const double pval = p[ij];
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double xp = pxval - pxsub;
            const double yp = pyval - pysub;
            const double zp = pzval - pzsub;
            const double pscale = scalep[ij];
            for (uint32_t kl = 0; kl < nkl; kl += 1) {
                const double qval = q[kl];
                const double pqmult = pval * qval;
                const double pqplus = pval + qval;
                const double pqpinv = 1. / pqplus;
                const double pqx = pxval - qx[kl];
                const double pqy = pyval - qy[kl];
                const double pqz = pzval - qz[kl];
                const double t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
                const double scale = pscale * scaleq[kl] / (pqmult * __builtin_sqrt(pqplus));
                struct Boys01 boys = boys01(t, scale);
                boys.f1 = -boys.f1 * qval * pqpinv;
                cbatch[0] += xp * boys.f0 + pqx * boys.f1;
                cbatch[1] += yp * boys.f0 + pqy * boys.f1;
                cbatch[2] += zp * boys.f0 + pqz * boys.f1;
            }
        }
    }
}
