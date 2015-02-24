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
         SUBROUTINE  OED__NAI_DERV_INT1D_TO_AB
     +
     +                    ( SHELLA,SHELLB,
     +                      NEXP,NGQPCEN,NGEXCEN,
     +                      NXYZA,NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_INT1D_TO_AB
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative nuclear attraction integrals [A|B], adding
C                up the contributions from all the 1D AB integrals
C                for all nuclear centers.
C
C                The routine uses the reduced Rys multiplication scheme
C                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889.
C                This scheme reuses intermediate products between
C                1DX and 1DY integrals.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x = A and B
C                    NEXP        =  current # of exponent pairs
C                    NGQPCEN     =  product of # of gaussian quadrature
C                                   points (roots) times # of nuclear
C                                   attraction centers
C                    NGEXCEN     =  product of # of gaussian quadrature
C                                   points times # of exponent pairs
C                                   times # of nuclear attraction
C                                   centers
C                    NXYZx       =  # of cartesian monomials for
C                                   x = A and B shells
C                    INT1Dx      =  all current 1D AB derivative
C                                   nuclear attraction integrals for
C                                   each cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   1D AB derivative integral products
C                    SCALE       =  the NGEXCEN scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian [A|B]
C                                   derivative nuclear attraction
C                                   integrals corresponding to all
C                                   current exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     DIFFY,DIFFZ

         INTEGER     I,J,K,M,N
         INTEGER     NEXP,NGQPCEN,NGEXCEN
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB
         INTEGER     XA,YA,ZA,XB,YB,ZB
         INTEGER     XAP,XBP
         INTEGER     YAMAX,YBMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGEXCEN)
         DOUBLE PRECISION  TEMP1 (1:NGEXCEN)
         DOUBLE PRECISION  TEMP2 (1:NGEXCEN)

         DOUBLE PRECISION  BATCH (1:NEXP,1:NXYZA,1:NXYZB)

         DOUBLE PRECISION  INT1DX (1:NGEXCEN,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INT1DY (1:NGEXCEN,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INT1DZ (1:NGEXCEN,0:SHELLA,0:SHELLB)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...outer loops over x,x-pairs. No skipping of
C                x,x-contribution of 0,0-type can be done here,
C                since the 1DX integrals carry the Rys weight!
C
C
         XBP = 0
         DO 100  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 110 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO N = 1,NGEXCEN
                  TEMP1 (N) = SCALE (N) * INT1DX (N,XA,XB)
               END DO
C
C
C             ...inner loops over y,y-pairs. Skip the multiplication
C                of y,y-contributions, if no y-coordinate derivative
C                was formed and we have a 0,0-pair, as then the 1DY
C                integrals are equal to 1.
C
C
               J = XBP
               DO 120 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 130 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGEXCEN
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGEXCEN
                            TEMP2 (N) = TEMP1 (N) * INT1DY (N,YA,YB)
                         END DO
                     END IF
C
C
C             ...skip multiplication of z,z-contributions, if we
C                have a 0,0-pair and no derivations were performed
C                on the z-coordinate, as then the 1DZ integrals
C                are equal to 1. All info concerning all three
C                x,y and z contributions have been collected for all
C                exponent pairs, nuclear centers and quadrature points
C                at once. Sum up the 1D X,Y,Z integral products over
C                the nuclear centers and the quadrature points to
C                the appropriate place of the [A|B] batch.
C
C
                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         K = 0
                         DO M = 1,NEXP
                            SUM = ZERO
                            DO N = 1,NGQPCEN
                               SUM = SUM + TEMP2 (K+N)
                            END DO
                            K = K + NGQPCEN
                            BATCH (M,I,J) = SUM
                         END DO
                     ELSE
                         K = 0
                         DO M = 1,NEXP
                            SUM = ZERO
                            DO N = 1,NGQPCEN
                               SUM = SUM + TEMP2 (K+N)
     +                                   * INT1DZ (K+N,ZA,ZB)
                            END DO
                            K = K + NGQPCEN
                            BATCH (M,I,J) = SUM
                         END DO
                     END IF
C
C
C             ...next z,z-pair and y,y-pair.
C
C
  130             CONTINUE
  120          CONTINUE
C
C
C             ...next x,x-pair.
C
C
               XAP = I
  110       CONTINUE
            XBP = J
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
