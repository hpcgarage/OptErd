C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU Lesser General Public License as published
C  by the Free Software Foundation; either version 2.1 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU Lesser General Public License for more details.

C  The GNU Lesser General Public License is included in this distribution
C  in the file COPYING.
         SUBROUTINE  OED__XYZ_E0_PCGTO_BLOCK
     +
     +                    ( NBATCH,NINT1D,
     +                      ATOMIC,
     +                      MIJ,NIJ,NIJBEG,NIJEND,
     +                      NPGTOA,NPGTOB,
     +                      NXYZET,NXYZP,
     +                      SHELLA,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      MOMENTX,MOMENTY,MOMENTZ,
     +                      ALPHAA,ALPHAB,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      XP,YP,ZP,
     +                      PAX,PAY,PAZ,PINVHF,SCALE,
     +                      INT1DX,SCRMTX,XMINUS,
     +                      INT1DY,SCRMTY,YMINUS,
     +                      INT1DZ,SCRMTZ,ZMINUS,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_E0_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__XYZ_1D_INTEGRALS
C                OED__XYZ_INT1D_TO_E0
C  DESCRIPTION : This operation calculates a batch of unnormed moment
C                integrals between primitive cartesian gaussians for
C                the shell range:
C
C                              [E|0]   , E = A to P = A + B
C                                   ij
C
C                and the block of ij exponent pairs. The total number
C                of moment integrals generated here is thus given by
C                the total number of cartesian monomials NXYZET times
C                the total number of exponent pairs MIJ in the present
C                block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,nxyzet)
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    integral batch
C                    NINT1D       =  space needed for each of the 1D
C                                    X,Y,Z integral arrays
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    NXYZET       =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,P=A+B
C                    NXYZP        =  # of cartesian monomials for
C                                    the P=A+B shell
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A and P=A+B
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B,P
C                    ABm          =  the m=x,y,z-coordinate differences
C                                    between centers A and B
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
C                    PRIMx        =  i,j labels of primitives for the
C                                    respective contraction shells
C                                    x = A,B
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B
C                    RHOAB        =  the complete set of NIJ exponential
C                                    prefactors between contraction
C                                    shells A and B
C                    PAx          =  will hold current MIJ coordinate
C                                    x=X,Y,Z differences P-A between
C                                    centers P and A
C                    PINVHF       =  will hold current MIJ values of
C                                    1/(2*P), where P are the exponent
C                                    sums for contraction shells A
C                                    and B
C                    SCALE        =  will hold current MIJ values of
C                                    scaling factors
C                    INT1Dx       =  will hold all current 1D integrals
C                                    for each cartesian component
C                                    (x = X,Y,Z)
C                    SCRMTx       =  scratch matrix to hold x-1 integrals
C                                    (x = X,Y,Z)
C                    xMINUS       =  scratch matrix that holds the (0|x-2|0)
C                                    integrals
C                    MOMENTx      =  tells which power of X,Y, or Z integrals
C                                    to compute
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian moment moment [E|0] 
C                                    integrals
C
C
C  AUTHOR      : Norbert Flocke
C                  - Wrote original OED package
C
C  MODIFIED    : Thomas Watson Jr.                    p  q  r
C                  - Modified OED package to compute X, Y, Z integrals
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC

         INTEGER     I,J,M
         INTEGER     IJ
         INTEGER     MIJ
         INTEGER     NBATCH,NINT1D
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NXYZET,NXYZP
         INTEGER     SHELLA,SHELLB,SHELLP
         INTEGER     MOMENTX,MOMENTY,MOMENTZ

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  PINV,PVAL
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  ZERO,HALF,ONE,ONEP5

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  XP (1:NIJ)
         DOUBLE PRECISION  YP (1:NIJ)
         DOUBLE PRECISION  ZP (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  XMINUS  (1:MIJ)
         DOUBLE PRECISION  YMINUS  (1:MIJ)
         DOUBLE PRECISION  ZMINUS  (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  SCALE   (1:MIJ)

         DOUBLE PRECISION  INT1DX  (1:NINT1D)
         DOUBLE PRECISION  INT1DY  (1:NINT1D)
         DOUBLE PRECISION  INT1DZ  (1:NINT1D)

         DOUBLE PRECISION  SCRMTX  (1:NINT1D)
         DOUBLE PRECISION  SCRMTY  (1:NINT1D)
         DOUBLE PRECISION  SCRMTZ  (1:NINT1D)

         PARAMETER  (ZERO  = 0.D0)
         PARAMETER  (HALF  = 0.5D0)
         PARAMETER  (ONE   = 1.D0)
         PARAMETER  (ONEP5 = 1.5D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the quantities needed to establish the
C                1D overlap integrals.
C
C
         IF (ATOMIC) THEN
             M = 0
             DO 100 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                PINV = ONE / (ALPHAA (I) + ALPHAB (J))
                XP (M) = XA
                YP (M) = YA
                ZP (M) = ZA
                PAX (M) = ZERO
                PAY (M) = ZERO
                PAZ (M) = ZERO
                PINVHF (M) = HALF * PINV
                SCALE  (M) = (PINV ** ONEP5) * NORMA (I) * NORMB (J)
  100        CONTINUE
         ELSE
             M = 0
             DO 110 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                EXPA = ALPHAA (I)
                EXPB = ALPHAB (J)
                PINV = ONE / (EXPA + EXPB)
                XP (M) = PINV * (EXPA * XA  +  EXPB * XB)
                YP (M) = PINV * (EXPA * YA  +  EXPB * YB)
                ZP (M) = PINV * (EXPA * ZA  +  EXPB * ZB)
                PVAL = - EXPB * PINV
                PAX (M) = PVAL * ABX
                PAY (M) = PVAL * ABY
                PAZ (M) = PVAL * ABZ
                PINVHF (M) = HALF * PINV
                SCALE  (M) = (PINV ** ONEP5)
     +                       * NORMA (I) * NORMB (J) * RHOAB (IJ)
  110        CONTINUE
         END IF
C
C
C             ...perform the following steps:
C
C                1) construct all 1D x,y,z moment integrals for
C                   all ij pairs.
C
C                2) assemble the complete [E|0] moment batch for all
C                   ij pairs using the 1D integrals. Arrays PAX and PAY
C                   are passed as scratch arrays.
C
C
                CALL    OED__XYZ_1D_INTEGRALS
     +
     +                    ( SHELLP,0,
     +                      ATOMIC,
     +                      MIJ,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      PAX,PAY,PAZ,
     +                      XP,YP,ZP,
     +                      PINVHF,
     +                      ZERO,ZERO,ZERO,
     +                      MOMENTX,MOMENTY,MOMENTZ,
     +                      XMINUS,YMINUS,ZMINUS,
     +                      SCRMTX,SCRMTY,SCRMTZ,
     +
     +                                INT1DX,
     +                                INT1DY,
     +                                INT1DZ )
     +
             CALL    OED__XYZ_INT1D_TO_E0
     +
     +                    ( SHELLA,SHELLP,
     +                      ATOMIC,
     +                      MIJ,
     +                      NXYZET,NXYZP,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      PAX,PAY,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
