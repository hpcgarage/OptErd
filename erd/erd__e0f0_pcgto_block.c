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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD_E0F0_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD_RYS_ROOTS_WEIGHTS */
/*                ERD_2D_COEFFICIENTS */
/*                ERD_2D_PQ_INTEGRALS */
/*                ERD_INT2D_TO_E000 */
/*                ERD_INT2D_TO_E0F0 */
/*                ERD_2D_ATOM_COEFFICIENTS */
/*                ERD_2D_ATOM_PQ_INTEGRALS */
/*                ERD_ATOM_INT2D_TO_E000 */
/*                ERD_ATOM_INT2D_TO_E0F0 */
/*  DESCRIPTION : This operation calculates a batch of unnormed electron */
/*                repulsion integrals between primitive cartesian */
/*                gaussians for the shell quadruplet range: */
/*                    [E0|F0]     , E = A to P, F = C to Q */
/*                           ijkl */
/*                and the block of ij and kl exponent pairs. The total */
/*                number of eris generated here is thus given by the */
/*                total number of cartesian monomials NXYZET*NXYZFT */
/*                times the total number of exponent pairs MIJKL in the */
/*                present block. */
/*                On exit, the batch elements will be stored as: */
/*                             batch (kl,ij,nxyzt) */
/*                  Input: */
/*                    NBATCH       =  size of the primitive cartesian */
/*                                    integral batch */
/*                    NINT2D       =  space needed for each of the 2D */
/*                                    X,Y,Z integral arrays */
/*                    ATOMIC       =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                    ATOMAB(CD)   =  indicates, if centers A and B */
/*                                    (C and D) coincide */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the contracted shell pairs A,B */
/*                                    (C,D) */
/*                    MIJKL        =  current # of ijkl primitive */
/*                                    index quadruplets (= MIJ*MKL) */
/*                    NIJ          =  total # of ij primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair A,B */
/*                    NIJBEG(END)  =  first(last) ij primitive index */
/*                                    defining the ij block */
/*                    NKL          =  total # of kl primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair C,D */
/*                    NKLBEG(END)  =  first(last) kl primitive index */
/*                                    defining the kl block */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    NGQSCR       =  size of gaussian quadrature */
/*                                    scratch space needed to calculate */
/*                                    all the quadrature roots */
/*                    MGQIJKL      =  # of roots times # of ijkl */
/*                                    quadruplets (= NGQP*MIJKL) */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = A,B,C,D */
/*                    NXYZE(F)T    =  sum of # of cartesian monomials */
/*                                    for all shells in the range */
/*                                    E = A,...,P=A+B and in the range */
/*                                    F = C,...,Q=C+D */
/*                    NXYZP(Q)     =  # of cartesian monomials for */
/*                                    the P=A+B and Q=C+D shells */
/*                    SHELLx       =  the shell type for contraction */
/*                                    shells x = A,P=A+B,C,Q=C+D */
/*                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers */
/*                                    x = A,B,C,D */
/*                    ABm(CDm)     =  the m=x,y,z-coordinate differences */
/*                                    between centers A and B (C and D) */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = A,B,C,D */
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
/*                                    x = A,B,C,D */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = A,B,C,D */
/*                    RHOAB(CD)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells A and B */
/*                                    (C and D) */
/*                    P            =  will hold current MIJ exponent */
/*                                    sums for contraction shells A */
/*                                    and B */
/*                    Px           =  will hold current MIJ coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers P=A+B */
/*                    PAx          =  will hold current MIJ coordinate */
/*                                    x=X,Y,Z differences P-A between */
/*                                    centers P and A */
/*                    PINVHF       =  will hold current MIJ values of */
/*                                    1/(2*P), where P are the exponent */
/*                                    sums for contraction shells A */
/*                                    and B */
/*                    SCALEP       =  will hold current MIJ values of */
/*                                    scaling factors related to point P */
/*                    Q            =  will hold current MKL exponent */
/*                                    sums for contraction shells C */
/*                                    and D */
/*                    Qx           =  will hold current MKL coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers Q=C+D */
/*                    QCx          =  will hold current MKL coordinate */
/*                                    x=X,Y,Z differences Q-C between */
/*                                    centers Q and C */
/*                    QINVHF       =  will hold current MKL values of */
/*                                    1/(2*Q), where Q are the exponent */
/*                                    sums for contraction shells C */
/*                                    and D */
/*                    SCALEQ       =  will hold current MKL values of */
/*                                    scaling factors related to point Q */
/*                    RTS          =  will hold all current MGQIJKL */
/*                                    quadrature roots */
/*                    WTS          =  will hold all current MGQIJKL */
/*                                    quadrature weights */
/*                    GQSCR        =  will be used as scratch space */
/*                                    for determining the quadrature */
/*                                    roots and weights */
/*                    TVAL         =  will hold current MIJKL values */
/*                                    of T-exponents defining the Rys */
/*                                    weight functions */
/*                    PQPINV       =  will hold current MIJKL values */
/*                                    of 1/(P+Q), i.e. the inverses */
/*                                    of all total exponent sums */
/*                    SCALEPQ      =  will hold current distinct MIJKL */
/*                                    (expanded to MGQIJKL) values of */
/*                                    the overal scaling factors for */
/*                                    the integrals */
/*                    Bxx          =  will hold the current MGQIJKL */
/*                                    coordinate independent VRR */
/*                                    B-coefficients (xx=00,01,10) */
/*                    C00x         =  will hold the current MGQIJKL */
/*                                    VRR C-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center P */
/*                    D00x         =  will hold the current MGQIJKL */
/*                                    VRR D-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center Q */
/*                    INT2Dx       =  will hold all current 2D integrals */
/*                                    for each cartesian component */
/*                                    (x = X,Y,Z) */
/*                  Output: */
/*                    BATCH        =  current batch of primitive */
/*                                    cartesian [E0|F0] integrals */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__e0f0_pcgto_block(
    uint32_t A, uint32_t B, uint32_t C, uint32_t D,
    uint32_t nij, uint32_t nkl,
    uint32_t nxyzet, uint32_t nxyzft,
    uint32_t nxyzp, uint32_t nxyzq,
    const uint32_t shell[restrict static 1],
    const double xyz0[restrict static 1],
    const double *restrict alpha[restrict static 1],
    const double *restrict cc[restrict static 1],
    int **vrrtab,
    const uint32_t prima[restrict static nij], const uint32_t primb[restrict static nij], const uint32_t primc[restrict static nkl], const uint32_t primd[restrict static nkl],
    const double norma[restrict static nij], const double normb[restrict static nij], const double normc[restrict static nkl], const double normd[restrict static nkl],
    const double rhoab[restrict static nij], const double rhocd[restrict static nkl],
    double output_buffer[restrict])
{
#ifdef __ERD_PROFILE__   
    #ifdef _OPENMP
    const int tid = omp_get_thread_num();
    #else
    const int tid = 0;
    #endif
#endif
    double xa = xyz0[A*4], ya = xyz0[A*4+1], za = xyz0[A*4+2];
    double xb = xyz0[B*4], yb = xyz0[B*4+1], zb = xyz0[B*4+2];
    double xc = xyz0[C*4], yc = xyz0[C*4+1], zc = xyz0[C*4+2];
    double xd = xyz0[D*4], yd = xyz0[D*4+1], zd = xyz0[D*4+2];
    const uint32_t shella = shell[A], shellb = shell[B], shellc = shell[C], shelld = shell[D];
    const uint32_t shellp = shella + shellb;
    const uint32_t shellq = shellc + shelld;
    const uint32_t shellt = shellp + shellq;
    const uint32_t ngqp = shellt / 2 + 1;
    
/*            ...predetermine 2D integral case. This is done in */
/*               order to distinguish the P- and Q-shell combinations */
/*               for efficient evaluation of the 2D integrals and */
/*               their VRR coefficients. The cases distinguished */
/*               are summarized in the following table, indicating */
/*               the value of CASE2D: */
/*                                Q-shell */
/*                              s    p   >p */
/*                            --------------- */
/*                           | */
/*                         s |  1    4    7 */
/*                           | */
/*                P-shell  p |  2    5    8 */
/*                           | */
/*                        >p |  3    6    9 */
/*                           | */   
/*             ...predetermine in 'K2' loops the quantities associated */
/*                with the A,B-part and C,D-part. Set the atom equality */
/*                case CASEAT here to exploit simplifications due to */
/*                center equalities later on: */
/*                     CASEAT = 1  -->    atomic (AA|AA) integrals */
/*                            = 2  -->  2-center (AA|CC) integrals */
/*                            = 3  -->  3-center (AB|CC) integrals */
/*                            = 4  -->  3-center (AA|CD) integrals */
/*                            = 5  -->  4-center (AB|CD) integrals */
    const uint32_t case2d = min32u(2, shellq) * 3 + min32u(2, shellp) + 1;
    const uint32_t nijkl = nij * nkl;
    const uint32_t mgqijkl = ngqp * nijkl;
    
    const double *restrict alphaa = alpha[A], *restrict alphab = alpha[B];
    const double *restrict cca = cc[A], *restrict ccb = cc[B];
    const size_t simd_nij = PAD_LEN(nij);
    ERD_SIMD_ALIGN double p[simd_nij], px[simd_nij], py[simd_nij], pz[simd_nij], pinvhf[simd_nij], scalep[simd_nij];
    ERD_SIMD_ZERO_TAIL_64f(p, simd_nij);
    ERD_SIMD_ZERO_TAIL_64f(px, simd_nij);
    ERD_SIMD_ZERO_TAIL_64f(py, simd_nij);
    ERD_SIMD_ZERO_TAIL_64f(pz, simd_nij);
    ERD_SIMD_ZERO_TAIL_64f(pinvhf, simd_nij);
    ERD_SIMD_ZERO_TAIL_64f(scalep, simd_nij);
    #pragma simd
    #pragma vector aligned
    for (uint32_t ij = 0; ij < nij; ++ij) {        
        const uint32_t i = prima[ij];
        const uint32_t j = primb[ij];
        const double expa = alphaa[i];
        const double expb = alphab[j];
        const double pval = expa + expb;
        const double pinv = 1.0 / pval;
        p[ij] = pval;
        px[ij] = (expa * xa + expb * xb) * pinv;
        py[ij] = (expa * ya + expb * yb) * pinv;
        pz[ij] = (expa * za + expb * zb) * pinv;
        pinvhf[ij] = pinv * 0.5;
        scalep[ij] = norma[i] * normb[j] * rhoab[ij] * cca[i] * ccb[j];
    }

    const double *restrict alphac = alpha[C], *restrict alphad = alpha[D];
    const double *restrict ccc = cc[C], *restrict ccd = cc[D];
    const size_t simd_nkl = PAD_LEN(nkl);
    ERD_SIMD_ALIGN double q[simd_nkl], qx[simd_nkl], qy[simd_nkl], qz[simd_nkl], qinvhf[simd_nkl], scaleq[simd_nkl];
    ERD_SIMD_ZERO_TAIL_64f(q, simd_nkl);
    ERD_SIMD_ZERO_TAIL_64f(qx, simd_nkl);
    ERD_SIMD_ZERO_TAIL_64f(qy, simd_nkl);
    ERD_SIMD_ZERO_TAIL_64f(qz, simd_nkl);
    ERD_SIMD_ZERO_TAIL_64f(qinvhf, simd_nkl);
    ERD_SIMD_ZERO_TAIL_64f(scaleq, simd_nkl);
    #pragma simd
    #pragma vector aligned
    for (uint32_t kl = 0; kl < nkl; ++kl) {
        const uint32_t k = primc[kl];
        const uint32_t l = primd[kl];
        const double expc = alphac[k];
        const double expd = alphad[l];
        double qval = expc + expd;
        const double qinv = 1.0 / qval;
        q[kl] = qval;
        qx[kl] = (expc * xc + expd * xd) * qinv;
        qy[kl] = (expc * yc + expd * yd) * qinv;
        qz[kl] = (expc * zc + expd * zd) * qinv;
        qinvhf[kl] = qinv * 0.5;
        scaleq[kl] = normc[k] * normd[l] * rhocd[kl] * ccc[k] * ccd[l];
    }

/*             ...the 'K4' loop over all ij- and kl-exponent pairs */
/*                in present ij and kl block to calculate all T's */
/*                and scaling factors for the cases: */
/*                     CASEAT = 1  -->    atomic (AA|AA) integrals */
/*                            = 2  -->  2-center (AA|CC) integrals */
/*                            = 3  -->  3-center (AB|CC) integrals */
/*                            = 4  -->  3-center (AA|CD) integrals */
/*                            = 5  -->  4-center (AB|CD) integrals */
/*                4-center (AB|CD) integrals are checked first */
/*                (most common occurence in large systems). */
    const size_t simd_mgqijkl = PAD_LEN(mgqijkl);
    const size_t simd_nijkl = PAD_LEN(nijkl);
    const uint32_t nint2d = simd_mgqijkl * (shellp + 1) * (shellq + 1);
    ERD_SIMD_ALIGN double int2dx[nint2d];
    ERD_SIMD_ALIGN double tval[simd_nijkl], pqpinv[simd_nijkl];
    ERD_SIMD_ZERO_TAIL_64f(tval, simd_nijkl);
    ERD_SIMD_ZERO_TAIL_64f(pqpinv, simd_nijkl);
    uint32_t m = 0;
    for (uint32_t ij = 0; ij < nij; ++ij) {
        const double pval = p[ij];
        const double pxval = px[ij];
        const double pyval = py[ij];
        const double pzval = pz[ij];
        const double pscale = scalep[ij];
        for (uint32_t kl = 0; kl < nkl; ++kl) {
            const double qval = q[kl];
            const double pqmult = pval * qval;
            const double pqplus = pval + qval;
            const double invers = 1.0 / pqplus;
            const double pqx = pxval - qx[kl];
            const double pqy = pyval - qy[kl];
            const double pqz = pzval - qz[kl];
            tval[m] = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * invers;
            pqpinv[m] = invers;
            int2dx[m] = pscale * scaleq[kl] / (pqmult * sqrt(pqplus));
            m++;
        }
    }

/*             ...if necessary, expand the scaling array size from */
/*                MIJKL to MGQIJKL starting from the last elements. */
    if (ngqp > 1) {
        uint32_t n = mgqijkl;
        for (uint32_t m = nijkl; m >= 1; m--) {
            for (uint32_t i = 1; i <= ngqp; i++) {
                int2dx[n - i] = int2dx[m - 1];
            }
            n -= ngqp;
        }
    }
    memset(&int2dx[mgqijkl], 0, sizeof(double) * (simd_mgqijkl - mgqijkl));


/*             ...calculate all roots and weights. */
    ERD_SIMD_ALIGN double rts[simd_mgqijkl];
    ERD_SIMD_ZERO_TAIL_64f(rts, simd_mgqijkl);
    const uint32_t nmom = (ngqp << 1) - 1;
    ERD_PROFILE_START(erd__rys_roots_weights)
    erd__rys_roots_weights(nijkl, ngqp, nmom, tval, rts, int2dx);
    ERD_PROFILE_END(erd__rys_roots_weights)
/*             ...perform the following steps: */
/*                1) generate all VRR coefficients. */
/*                2) construct all 2D PQ x,y,z integrals using all the */
/*                   weights and all the generated VRR coefficients for */
/*                   all exponent quadruples. */
/*                3) assemble the complete [E0|F0] batch for all ij and */
/*                   kl pairs using the 2D integrals. Arrays B00 and B01 */
/*                   are passed as scratch arrays. */
/*                The last step 3) is the most compute intensive and */
/*                separate routines are provided depending on presence */
/*                of s-shells. The gainings are in the innermost loops */
/*                of these routines, which are considerably simplified */
/*                for the special s-shell cases. Note, that the case */
/*                in which both P- and Q-shells are s-shells cannot */
/*                arise, as this case is dealt with in separate routines. */
    ERD_SIMD_ALIGN double b00[simd_mgqijkl], b01[simd_mgqijkl], b10[simd_mgqijkl], c00x[simd_mgqijkl], c00y[simd_mgqijkl], c00z[simd_mgqijkl], d00x[simd_mgqijkl], d00y[simd_mgqijkl], d00z[simd_mgqijkl];
    ERD_SIMD_ZERO_TAIL_64f(b00, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(b01, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(b10, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(c00x, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(c00y, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(c00z, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(d00x, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(d00y, simd_mgqijkl);
    ERD_SIMD_ZERO_TAIL_64f(d00z, simd_mgqijkl);
    ERD_PROFILE_START(erd__2d_coefficients)
    erd__2d_coefficients(nij, nkl, ngqp, p, q,
                          px, py, pz, qx, qy, qz,
                          &xyz0[A*4], &xyz0[C*4],
                          pinvhf, qinvhf, pqpinv, rts,
                          b00, b01, b10,
                          c00x, c00y, c00z,
                          d00x, d00y, d00z);
    ERD_PROFILE_END(erd__2d_coefficients)

    ERD_SIMD_ALIGN double int2dy[nint2d], int2dz[nint2d];
    ERD_PROFILE_START(erd__2d_pq_integrals)
    erd__2d_pq_integrals(shellp, shellq, simd_mgqijkl,
                          b00, b01, b10, c00x, c00y, c00z, d00x,
                          d00y, d00z, case2d,
                          int2dx, int2dy, int2dz);
    ERD_PROFILE_END(erd__2d_pq_integrals)
    
    ERD_PROFILE_START(erd__int2d_to_e0f0)
    erd__int2d_to_e0f0(shella, shellp, shellc, shellq,
                        simd_mgqijkl, nxyzet, nxyzft,
                        int2dx, int2dy, int2dz, vrrtab, output_buffer);
    ERD_PROFILE_END(erd__int2d_to_e0f0)
}
