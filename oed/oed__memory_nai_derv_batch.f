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
         SUBROUTINE  OED__MEMORY_NAI_DERV_BATCH
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      NUCLEI,
     +                      XN,YN,ZN,
     +                      IXDERC,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      DERCX,DERCY,DERCZ,
     +                      ALPHA,CC,
     +                      SPHERIC,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__MEMORY_NAI_DERV_BATCH
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__MEMORY_NAI_DERV_CSGTO
C  DESCRIPTION : Main operation that calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted electron nuclear attraction derivative
C                integrals.
C
C
C                  Input (x = 1 and 2):
C
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                    NUCLEI       =  # of nuclear attraction centers
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    IXDERC       =  the index of which of the nuclear
C                                    attraction centers is to be
C                                    differentiated. If that index
C                                    corresponds to one of the centers
C                                    1 and/or 2, it will already be
C                                    differentiated along with these
C                                    centers, hence values transmitted
C                                    for DERCX,DERCY,DERCZ are
C                                    irrelevant in that case. If no
C                                    nuclear attraction center is to
C                                    be differentiated besides those
C                                    which are possibly equal to
C                                    centers 1 and/or 2, set this
C                                    index value =< 0
C                    DERyp        =  the order of differentiation on
C                                    centers y = 1 and 2 with respect
C                                    to the p = x,y,z coordinates
C                    DERCp        =  the order of differentiation for
C                                    the IXDERC-th nuclear attraction
C                                    center with respect to the
C                                    p = x,y,z coordinates
C                    ALPHA        =  primitive exponents for csh
C                                    1 and 2 in that order
C                    CC           =  contraction coefficient for csh
C                                    1 and 2 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C
C                  Output:
C
C                    IMIN,IOPT    =  minimum/optimum integer memory
C                    ZMIN,ZOPT    =  minimum/optimum flp memory
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INCLUDE     'oed__tuning.inc'

         LOGICAL     SPHERIC

         INTEGER     DER1X,DER1Y,DER1Z
         INTEGER     DER2X,DER2Y,DER2Z
         INTEGER     DERCX,DERCY,DERCZ
         INTEGER     IMIN,IOPT
         INTEGER     IXDERC
         INTEGER     NALPHA,NCOEFF
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NUCLEI
         INTEGER     SHELL1,SHELL2
         INTEGER     ZMIN,ZOPT

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2

         DOUBLE PRECISION  XN (1:NUCLEI)
         DOUBLE PRECISION  YN (1:NUCLEI)
         DOUBLE PRECISION  ZN (1:NUCLEI)

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
C
C
C------------------------------------------------------------------------
C
C
C             ...call memory routine.
C
C
         CALL  OED__MEMORY_NAI_DERV_CSGTO
     +
     +              ( NALPHA,NCOEFF,
     +                NCGTO1,NCGTO2,
     +                NPGTO1,NPGTO2,
     +                SHELL1,SHELL2,
     +                X1,Y1,Z1,X2,Y2,Z2,
     +                NUCLEI,
     +                XN,YN,ZN,
     +                IXDERC,
     +                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DERCX,DERCY,DERCZ,
     +                ALPHA,CC,
     +                L1CACHE,NCTROW,
     +                SPHERIC,
     +
     +                         IMIN,IOPT,
     +                         ZMIN,ZOPT )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
