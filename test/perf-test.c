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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <sys/time.h>
#include <config.h>
#include <inttypes.h>

#include <flops.h>
#include <ratios.h>
#include <screening.h>
#include <erd_profile.h>

static uint64_t get_cpu_frequency(void) {
    const uint64_t start_clock = __rdtsc();
    sleep(1);
    const uint64_t end_clock = __rdtsc();
    return end_clock - start_clock;
}

static inline uint32_t xorshift_rand(uint32_t* state) {
    uint32_t y = *state;
    y ^= y << 13;
    y ^= y >> 17;
    y ^= y << 5;
    *state = y;
    return y;
}

int main (int argc, char **argv)
{
    int nnz;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    if (argc != 5) {
        printf ("Usage: %s <basisset> <xyz> <fraction> <nthreads>\n", argv[0]);
        return -1;
    }

    const uint64_t freq = get_cpu_frequency();

    const double fraction = atof(argv[3]);
    assert(fraction > 0.0 && fraction <= 1.0);
    const int nthreads = atoi(argv[4]);
    #ifdef _OPENMP
    omp_set_num_threads(nthreads);
    #else
    assert(nthreads == 1);
    #endif
    
    // load basis set
    BasisSet_t basis;
    CInt_createBasisSet(&basis);
    CInt_loadBasisSet(basis, argv[1], argv[2]);
    schwartz_screening(basis, &shellptr, &shellid, &shellrid, &shellvalue, &nnz);

    printf("Molecule info:\n");
    printf("  #Atoms\t= %d\n", CInt_getNumAtoms(basis));
    printf("  #Shells\t= %d\n", CInt_getNumShells(basis));
    printf("  #Funcs\t= %d\n", CInt_getNumFuncs(basis));
    printf("  #OccOrb\t= %d\n", CInt_getNumOccOrb(basis));
    printf("  nthreads\t= %d\n", nthreads);

    ERD_t erd;
    CInt_createERD(basis, &erd, nthreads);

    double* totalcalls = (double *) malloc(sizeof (double) * nthreads * 64);
    assert(totalcalls != NULL);
    double* totalnintls = (double *) malloc(sizeof (double) * nthreads * 64);
    assert(totalnintls != NULL);

    #pragma omp parallel for
    for (int i = 0; i < nthreads; i++) {
        totalcalls[i * 64] = 0.0;
        totalnintls[i * 64] = 0.0;
    }
    printf("Computing integrals ...\n");

    // reset profiler
    erd_reset_profile();

    //printf ("max memory footprint per thread = %lf KB\n",
    //    CInt_getMaxMemory (erd[0])/1024.0);

    const uint32_t shellCount = CInt_getNumShells(basis);

    /* In (fraction) cases rand() returns value not greater than computationThreshold */
    const uint32_t computationThresholdMax = 0xFFFFFFFEu;
    const uint32_t computationThreshold = lround(fraction * computationThresholdMax);

    const uint64_t start_clock = __rdtsc();
        
    #pragma omp parallel
    {
        #ifdef _OPENMP
        const int tid = omp_get_thread_num();
        #else
        const int tid = 0;
        #endif

        BEGIN_RECORD_FLOPS
        BEGIN_RECORD_RATIO

        #pragma omp for schedule(dynamic)
        for (uint32_t shellIndexMP = 0; shellIndexMP < shellCount * shellCount; shellIndexMP++) {
            const uint32_t shellIndexM = shellIndexMP / shellCount;
            const uint32_t shellIndexP = shellIndexMP % shellCount;
            
            const uint32_t shellIndexNStart = shellptr[shellIndexM];
            const uint32_t shellIndexNEnd = shellptr[shellIndexM+1];
            const uint32_t shellIndexQStart = shellptr[shellIndexP];
            const uint32_t shellIndexQEnd = shellptr[shellIndexP+1];

            /* Should depend only on loop iteration to guarantee the same value in optimized and reference versions regardless of thread order */
            uint32_t rng_state = shellIndexMP;
            /* Do several iteration to make starting value less uniform */
            xorshift_rand(&rng_state);
            xorshift_rand(&rng_state);
            xorshift_rand(&rng_state);

            /* Prepare indices */
            for (uint32_t shellIndexNOffset = shellIndexNStart; shellIndexNOffset != shellIndexNEnd; shellIndexNOffset++) {
                const uint32_t shellIndexN = shellid[shellIndexNOffset];
                if (shellIndexM > shellIndexN)
                    continue;

                const double shellValueMN = shellvalue[shellIndexNOffset];
                for (uint32_t shellIndexQOffset = shellIndexQStart; shellIndexQOffset != shellIndexQEnd; shellIndexQOffset++) {
                    const uint32_t shellIndexQ = shellid[shellIndexQOffset];
                    if (shellIndexP > shellIndexQ)
                        continue;

                    if (shellIndexM + shellIndexN > shellIndexP + shellIndexQ)
                        continue;

                    /* Sample random integer. With probability (fraction) process the shell quartet. */
                    if ((computationThreshold == computationThresholdMax) || (xorshift_rand(&rng_state) <= computationThreshold)) {
                        double *integrals;
                        int nints;
                        CInt_computeShellQuartet(basis, erd, tid, shellIndexM, shellIndexN, shellIndexP, shellIndexQ, &integrals, &nints);

                        totalcalls[tid * 64] += 1;
                        totalnintls[tid * 64] += nints;
                    }
                }
            }
        }
        
        END_RECORD_FLOPS
        END_RECORD_RATIO
    }
    const uint64_t end_clock = __rdtsc();

    for (int i = 1; i < nthreads; i++) {
        totalcalls[0 * 64] = totalcalls[0 * 64] + totalcalls[i * 64];
        totalnintls[0 * 64] = totalnintls[0 * 64] + totalnintls[i * 64];
    } 
    const uint64_t total_ticks = end_clock - start_clock;
    const double timepass = ((double) total_ticks) / freq;

    printf("Done\n");
    printf("\n");
    printf("Number of calls: %.6le, Number of integrals: %.6le\n", totalcalls[0], totalnintls[0]);
    printf("Total GigaTicks: %.3lf, freq = %.3lf GHz\n", (double) (total_ticks) * 1.0e-9, (double)freq/1.0e9);
    printf("Total time: %.4lf secs\n", timepass);
    printf("Average time per call: %.3le us\n", 1000.0 * 1000.0 * timepass / totalcalls[0]);

    // use 1 if thread timing is not required
    erd_print_profile(1);

    REPORT_FLOPS
    REPORT_RATIO

    CInt_destroyERD(erd);
    free(totalcalls);
    free(totalnintls);
    free(shellptr);
    free(shellid);
    free(shellvalue);
    free(shellrid);

    CInt_destroyBasisSet(basis);

    return 0;
}
