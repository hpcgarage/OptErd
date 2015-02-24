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
         SUBROUTINE  OED__NAI_E0_DEF_BLOCKS
     +
     +               ( ZMAX,
     +                 NPGTOA,NPGTOB,
     +                 SHELLP,
     +                 NIJ,NRS,NCEN,
     +                 NGQP,NGQSCR,
     +                 NXYZT,
     +                 L1CACHE,NCTROW,
     +                 MEMORY,
     +
     +                     NIJBLK,
     +                     NPSIZE,NCSIZE,NWSIZE,
     +                     NINT1D,
     +                     MXPRIM,MNPRIM,
     +                     ZCBATCH,ZPBATCH,ZWORK,
     +                     ZNORMA,ZNORMB,
     +                     ZRHOAB,
     +                     ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                     ZRTS,ZWTS,ZGQSCR,ZTVAL,
     +                     ZR1X,ZR1Y,ZR1Z,ZR2,
     +                     ZINT1DX,ZINT1DY,ZINT1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_E0_DEF_BLOCKS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted (e|0) nuclear
C                attraction batch. It determines the block size
C                partitions (if any) for the ij exponent paris of the
C                whole primitive [e|0] batch and returns pointers to
C                the flp data sections needed by the (e|0) generation.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (e|0) generation flp
C                pointers. This is useful for establishing just the
C                overall memory size needed for the (e|0) generation.
C                The keyword MEMORY has to be set true for this case
C                and obviously the value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=A,B
C                    SHELLP       =  the shell sum A+B
C                    NIJ          =  total # of ij primitive index
C                                    pairs corresponding to the csh
C                                    pairs A,B
C                    NRS          =  total # of rs contraction index
C                                    pairs corresponding to the csh
C                                    pairs A,B
C                    NCEN         =  # of nuclear attraction centers
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NGQSCR       =  size of gaussian quadrature
C                                    scratch space needed to calculate
C                                    all the quadrature roots
C                    NXYZT        =  total # of cartesian monomial
C                                    pairs
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    MEMORY       =  if this keyword is true, the
C                                    routine will only determine the
C                                    minimum / optimum flp memory and
C                                    store these values into NPSIZE and
C                                    NCSIZE, respectively (see below)
C
C                  Output:
C
C                    NIJBLK       =  contains the block size for the
C                                    ij primitive indices in order to
C                                    perform efficient contractions.
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays.
C                                    If MEMORY is true, NPSIZE and
C                                    NCSIZE will contain respectively
C                                    the minimum and optimum flp memory
C                    NINT1D       =  space needed for the 1D X,Y,Z
C                                    integral arrays
C                    MXPRIM       =  the maximum # of primitives
C                                    between all i and j primitives,
C                                    i.e. = max (i,j)
C                    MNPRIM       =  the minimum # of primitives
C                                    between all i and j primitives,
C                                    i.e. = min (i,j)
C                    Z.....       =  pointers for space partition of
C                                    the flp array (see below)
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  |
C
C
C                   Zone 1: final (e|0) batch
C                   Zone 2: block [e|0]
C                   Zone 3 : complete set of A,B norms and complete
C                            set (after screening!) of exponential
C                            prefactors and nuclear attraction center
C                            information (coordinates + nuclear charges)
C                   Zone 4 : i) scratch for block [e|0] generation
C                           ii) scratch for (e|0) generation
C
C
C                Memory allocation offsets for the primitive [e|0]
C                batches generation + contraction to (e|0):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (e|0) batch
C
C                  --- Zone 2 ---
C
C                   ZPBATCH = offset for (blocked) [e|0] batch
C
C                  --- Zone 3 ---
C
C                   ZNORMA = offset for all A norms
C                   ZNORMB = offset for all B norms
C                   ZRHOAB = offset for all AB exponential prefactors
C
C                  --- Zone 4: for block [e|0] generation only ---
C
C                   ZPX = offset for (blocked) PX values
C                   ZPY = offset for (blocked) PY values
C                   ZPZ = offset for (blocked) PZ values
C                   ZPAX = offset for (blocked) PAX values
C                   ZPAY = offset for (blocked) PAY values
C                   ZPAZ = offset for (blocked) PAZ values
C                   ZPINVHF = offset for (blocked) 1/2P values
C                   ZSCALE = offset for (blocked) scale values
C
C                   ZRTS = offset for (blocked) quadrature roots 
C                   ZWTS = offset for (blocked) quadrature weights
C                   ZGQSCR = offset for (blocked) quadrature scratch
C                   ZTVAL = offset for (blocked) T exponents
C
C                   ZR1X = offset for (blocked) VRR coeff R1X
C                   ZR1Y = offset for (blocked) VRR coeff R1Y
C                   ZR1Z = offset for (blocked) VRR coeff R1Z
C                   ZR2 = offset for (blocked) VRR coeff R2
C
C                   ZINT1DX = offset for (blocked) 1DX integrals
C                   ZINT1DY = offset for (blocked) 1DY integrals
C                   ZINT1DZ = offset for (blocked) 1DZ integrals
C
C                  --- Zone 4: for contraction only ---
C
C                   ZWORK = offset for contraction working array
C
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     MEMORY

         INTEGER     IJDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     N
         INTEGER     MIJ,MIJBIG,NCEN,MIJCEN,MGIJCEN
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NGQP,NGQSCR
         INTEGER     NIJ,NIJBLK
         INTEGER     NINT1D
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ
         INTEGER     SHELLP
         INTEGER     WCOL
         INTEGER     ZMAX
         INTEGER     ZONE12,ZONE3,ZONE4,ZONE4B,ZONE4C
         INTEGER     ZONEMIN,ZONEOPT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZRHOAB,
     +               ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZRTS,ZWTS,ZGQSCR,ZTVAL,
     +               ZR1X,ZR1Y,ZR1Z,ZR2,
     +               ZINT1DX,ZINT1DY,ZINT1DZ

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [e|0] blocks:
C
C                The next few lines of code are among the most
C                important for performance of the entire overlap
C                integral package. Any badly chosen sizes for the
C                primitive blocks has a severe effect on the integral
C                evaluation timing. Two opposite effects have
C                to be considered:
C
C                1) Reuse level 1 cache data as much as possible,
C                   i.e. choose the primitive blocks in such a
C                   way that most of the heavy computational work
C                   is being done reusing data in each cache line.
C
C                2) Choose the size of the primitive blocks as
C                   large as possible observing point 1), to reduce
C                   subroutine call overheads as much as possible.
C
C                After some experimenting, the following procedure
C                was found to work best:
C
C                  i) The optimum [e|0] block size is assumed to
C                     be directly proportional to the level 1 cache
C                     size.
C
C                 ii) The proportionality factor PFACT has to be
C                     determined by 'experiment', running a real
C                     molecular case with fully contracted basis sets.
C
C
         PFACT = 4

         MIJ = PFACT * L1CACHE / NXYZT
         MIJ = MAX (1,MIJ)
         MIJ = MIN (MIJ,NIJ)
C
C
C             ...the optimum block size MIJ has been determined.
C                We have to see now how it fits with the rest of
C                the needed data into the maximum memory given.
C                If necessary subdivide the optimum block size
C                into convenient subblocks that fit into memory.
C
C                The subdivision of the MIJ block is done in such
C                a way that successive powers of 2 are checked as
C                divisors of MIJ. The maximum # of division steps
C                MXSTEP can thus be calculated beforehand by knowing
C                the fact that the smallest subdivided MIJ block
C                is of size 1. The succesive divisors are checked
C                in the following order:
C
C
C                            step  |  div of MIJ
C                           ----------------------
C                              1   |       1
C                              2   |       2
C                              3   |       4
C                              4   |       8
C                              5   |      16
C                             ...         ...
C
C                The routine stops, if even a MIJ block of minimum
C                size 1 cannot be accomodated into the memory.
C
C
         MXPRIM = MAX (NPGTOA,NPGTOB)
         MNPRIM = MIN (NPGTOA,NPGTOB)

         NCSIZE  = NXYZT * NRS

         ZONE3 = NPGTOA + NPGTOB + NIJ
C
C
C             ...if the MEMORY keyword is activated, the routine
C                will only determine the optimum and minimum flp
C                memory (in that order), place them respectively
C                into the NCSIZE and NPSIZE variables and exit.
C
C
         IF (MEMORY) THEN

             MIJCEN  = MIJ * NCEN
             MGIJCEN = NGQP * MIJCEN
             NINT1D  = MGIJCEN * (SHELLP+1)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = NGQSCR + MIJCEN + 7*(MIJ+MGIJCEN) + 3*NINT1D
             ZONE4C  = NWSIZE
             ZONE4   = MAX (ZONE4B,ZONE4C)

             ZONEOPT = ZONE12 + ZONE3 + ZONE4

             MIJ     = 1
             MIJCEN  = MIJ * NCEN
             MGIJCEN = NGQP * MIJCEN
             NINT1D  = MGIJCEN * (SHELLP+1)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = NGQSCR + MIJCEN + 7*(MIJ+MGIJCEN) + 3*NINT1D
             ZONE4C  = NWSIZE
             ZONE4   = MAX (ZONE4B,ZONE4C)

             ZONEMIN = ZONE12 + ZONE3 + ZONE4

             NPSIZE  = ZONEMIN
             NCSIZE  = ZONEOPT

             RETURN

         ELSE
C
C
C             ...the actual fitting into the maximum memory given.
C
C
             POW2IJ = DLOG ( DFLOAT (MIJ) ) / DLOG2
             MXSTEP = POW2IJ + 1

             IJDIV = 1
             MIJBIG = MIJ

             DO 100 N = 1,MXSTEP

                MIJ = MAX (1, MIJBIG / IJDIV )

                MIJCEN  = MIJ * NCEN
                MGIJCEN = NGQP * MIJCEN
                NINT1D  = MGIJCEN * (SHELLP+1)
                NPSIZE  = NXYZT * MIJ
                WCOL    = MNPRIM
                NWSIZE  = NCTROW * WCOL
                NPSIZE  = MAX (NPSIZE,NWSIZE)
                NWSIZE  = NPSIZE
                ZONE12  = NPSIZE + NCSIZE
                ZONE4B  = NGQSCR + MIJCEN + 7*(MIJ+MGIJCEN) + 3*NINT1D
                ZONE4C  = NWSIZE
                ZONE4   = MAX (ZONE4B,ZONE4C)

                IF (ZONE12+ZONE3+ZONE4.LE.ZMAX) THEN

                    NIJBLK = MIJ
                    NWSIZE = ZMAX - ZONE12 - ZONE3
C
C
C             ...generate the memory allocation pointers.
C
C
                    ZCBATCH = 1
                    ZPBATCH = NCSIZE + 1

                    ZNORMA = ZPBATCH + NPSIZE
                    ZNORMB = ZNORMA + NPGTOA
                    ZRHOAB = ZNORMB + NPGTOB

                    ZPX = ZRHOAB + NIJ
                    ZPY = ZPX + MIJ
                    ZPZ = ZPY + MIJ
                    ZPAX = ZPZ + MIJ
                    ZPAY = ZPAX + MIJ
                    ZPAZ = ZPAY + MIJ
                    ZPINVHF = ZPAZ + MIJ
                    ZSCALE = ZPINVHF + MIJ

                    ZRTS = ZSCALE + MGIJCEN
                    ZWTS = ZRTS + MGIJCEN
                    ZGQSCR = ZWTS + MGIJCEN
                    ZTVAL = ZGQSCR + NGQSCR

                    ZR1X = ZTVAL + MIJCEN
                    ZR1Y = ZR1X + MGIJCEN
                    ZR1Z = ZR1Y + MGIJCEN
                    ZR2 = ZR1Z + MGIJCEN

                    ZINT1DX = ZR2 + MGIJCEN
                    ZINT1DY = ZINT1DX + NINT1D
                    ZINT1DZ = ZINT1DY + NINT1D

                    ZWORK = ZPX

                    RETURN
                END IF

                IJDIV = IJDIV + IJDIV

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for (e|0) ! '
             WRITE (*,*) ' NIJ,MIJ = ',NIJ,MIJ
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (oed__nai_e0_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         END
