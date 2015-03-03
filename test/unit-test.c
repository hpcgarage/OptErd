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
#include <inttypes.h>
#include <stdint.h>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <sys/time.h>
#else /* Linux, etc. */
#include <time.h>
#endif

#ifdef HAS_MALLOC_H
#include <malloc.h>
#endif

#include <screening.h>

static uint64_t getTimerTicks(void) {
#if defined(__APPLE__) && defined(__MACH__)
    struct timeval tv;
    const int result = gettimeofday(&tv, NULL);
    if (result != 0) {
        perror("gettimeofday");
        abort();
    }
    return ((uint64_t)tv.tv_sec) * UINT64_C(1000000) + ((uint64_t)tv.tv_usec);
#else /* Linux, etc. */
    struct timespec ts;
    const int result = clock_gettime(CLOCK_MONOTONIC, &ts);
    if (result != 0) {
        perror("clock_gettime");
        abort();
    }
    return ((uint64_t)ts.tv_sec) * UINT64_C(1000000000) + ((uint64_t)ts.tv_nsec);
#endif
}

static uint64_t getTimerFrequency(void) {
    return UINT64_C(1000000000);
}

int main(int argc, char **argv) {

    if (argc != 3) {
        printf("Usage: %s <basisset> <xyz>\n", argv[0]);
        return 0;
    }
    FILE *ref_data_file = fopen ("ivalues.ref", "r");
    if (ref_data_file == NULL) {
        fprintf(stderr, "ivalues.ref does not exist\n");
        exit (0);
    }

    BasisSet_t basis;
    CInt_createBasisSet(&basis);
    CInt_loadBasisSet(basis, argv[1], argv[2]);
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    int nnz;
    schwartz_screening(basis, &shellptr, &shellid, &shellrid, &shellvalue, &nnz);

    printf("Molecule info:\n");
    printf("  #Atoms\t= %d\n", CInt_getNumAtoms(basis));
    printf("  #Shells\t= %d\n", CInt_getNumShells(basis));
    printf("  #Funcs\t= %d\n", CInt_getNumFuncs(basis));
    printf("  #OccOrb\t= %d\n", CInt_getNumOccOrb(basis));

    ERD_t erd;
    CInt_createERD(basis, &erd, 1);

    printf("Computing integrals ...\n");
    const uint32_t shellCount = CInt_getNumShells(basis);
    int dimMax = CInt_getMaxShellDim(basis);
    double *referenceIntegrals = (double *)malloc(sizeof(double) * dimMax * dimMax * dimMax * dimMax);
    assert(referenceIntegrals != NULL);
#if 0
    int nshells = shellCount;
    for (int M = 0; M < nshells; M++)
    {
        for (int N = 0; N < nshells; N++)
        {
            for (int P = 0; P < nshells; P++)
            {
                for (int Q = 0; Q < nshells; Q++)
                {
                    if (M < N || P < Q || (M+N) < (P+Q))
                        continue;
                    double *integrals;
                    int nints;
                    CInt_computeShellQuartet (basis, erd, 0,
                                              M, N, P, Q, &integrals, &nints);
                    int dimM = CInt_getShellDim (basis, M);
                    int dimN = CInt_getShellDim (basis, N);
                    int dimP = CInt_getShellDim (basis, P);
                    int dimQ = CInt_getShellDim (basis, Q);
                    int startM = CInt_getFuncStartInd (basis, M);
                    int startN = CInt_getFuncStartInd (basis, N);
                    int startP = CInt_getFuncStartInd (basis, P);
                    int startQ = CInt_getFuncStartInd (basis, Q);
                    for (int iM = 0; iM < dimM; iM++)
                    {
                        for (int iN = 0; iN < dimN; iN++)
                        {
                            for (int iP = 0; iP < dimP; iP++)
                            {
                                for (int iQ = 0; iQ < dimQ; iQ++)
                                {
                                    double ivalue = integrals[iM + dimM * (iN + dimN * (iP + dimP * iQ))];
                                    printf ("%d %d %d %d %le\n",
                                            startM + iM, startN + iN,
                                            startP + iP, startQ + iQ,
                                            ivalue);
                                }
                            }
                        }
                   }
               }
            }
        }
    }
#endif
    uint64_t totalCalls = 0;
    uint64_t totalIntegralsCount = 0;
    uint64_t totalTicks = 0;
    uint32_t errcount = 0;
    for (uint32_t shellIndexM = 0; shellIndexM != shellCount; shellIndexM++) {
        const uint32_t shellIndexNStart = shellptr[shellIndexM];
        const uint32_t shellIndexNEnd = shellptr[shellIndexM+1];
        for (uint32_t shellIndexP = 0; shellIndexP != shellCount; shellIndexP++) {
            const uint32_t shellIndexQStart = shellptr[shellIndexP];
            const uint32_t shellIndexQEnd = shellptr[shellIndexP+1];

            uint32_t* shellIndicesN = memalign(64, sizeof(double) * (shellIndexNEnd - shellIndexNStart) * (shellIndexQEnd - shellIndexQStart));
            uint32_t* shellIndicesQ = memalign(64, sizeof(double) * (shellIndexNEnd - shellIndexNStart) * (shellIndexQEnd - shellIndexQStart));

            uint32_t shellIndicesCount = 0;

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

                    const double shellValuePQ = shellvalue[shellIndexQOffset];
                    if (shellValueMN * shellValuePQ < TOLSRC * TOLSRC)
                        continue;

                    shellIndicesN[shellIndicesCount] = shellIndexN;
                    shellIndicesQ[shellIndicesCount] = shellIndexQ;
                    shellIndicesCount++;
                }
            }

            /* Validate results */
            for (uint32_t shellIndexIndex = 0; shellIndexIndex != shellIndicesCount; shellIndexIndex++) {
                /* Process batch of indices */
                totalCalls += 1;
                const uint64_t startTicks = getTimerTicks();

                double *integrals;
                int integralsCount;
                CInt_computeShellQuartets(basis, erd, 0 /* Thread ID */,
                    shellIndexM, &shellIndicesN[shellIndexIndex], shellIndexP, &shellIndicesQ[shellIndexIndex], 1,
                    &integrals, &integralsCount);

                const uint64_t endTicks = getTimerTicks();
                totalTicks += endTicks - startTicks;
                totalIntegralsCount += integralsCount;

                const uint32_t shellIndexN = shellIndicesN[shellIndexIndex];
                const uint32_t shellIndexQ = shellIndicesQ[shellIndexIndex];
                
                int referenceIntegralsCount;
                fread(&referenceIntegralsCount, sizeof(int), 1, ref_data_file);
                if (referenceIntegralsCount != 0) {
                    fread(referenceIntegrals, sizeof(double), referenceIntegralsCount, ref_data_file);
                }
                
                for (int i = 0; i < integralsCount; i++) {
                    if (isnan(integrals[i])) {
                        printf("ERROR: NAN INTEGRAL: %"PRIu32" %"PRIu32" %"PRIu32" %"PRIu32"\n",
                                shellIndexM, shellIndexN, shellIndexP, shellIndexQ);
                        errcount++;
                    }
                }

                if (integralsCount == 0 && referenceIntegralsCount == 0) {
                    continue;
                } else if (integralsCount != 0 && referenceIntegralsCount == 0) {
                    for (int k = 0; k < integralsCount; k++) {
                        if (integrals[k] > 1.0e-10) {
                            printf("ERROR: %"PRIu32" %"PRIu32" %"PRIu32" %"PRIu32": %le %le\n",
                                shellIndexM, shellIndexN, shellIndexP, shellIndexQ, 0.0, integrals[k]);
                            errcount++;
                        }
                    }
                } else if (integralsCount == 0 && referenceIntegralsCount != 0) {
                    for (int k = 0; k < referenceIntegralsCount; k++) {
                        if (referenceIntegrals[k] > 1e-10) {
                            printf("ERROR: %"PRIu32" %"PRIu32" %"PRIu32" %"PRIu32": %le %le\n",
                                shellIndexM, shellIndexN, shellIndexP, shellIndexQ, referenceIntegrals[k], 0.0);
                            errcount++;
                        }
                    }
                } else if (integralsCount == referenceIntegralsCount && integralsCount != 0) {
                    for (int k = 0; k < referenceIntegralsCount; k++) {
                        if ((fabs(referenceIntegrals[k]) < 1.0e-6) ||
                            (fabs(integrals[k]) < 1.0e-6))
                        {
                            if (fabs(referenceIntegrals[k] - integrals[k]) > 1.0e-10) {
                                printf ("1 ERROR: %"PRIu32" %"PRIu32" %"PRIu32" %"PRIu32": %le %le\n",
                                    shellIndexM, shellIndexN, shellIndexP, shellIndexQ, referenceIntegrals[k], integrals[k]);
                                errcount++;
                            }
                        } else {
                            if ((fabs(referenceIntegrals[k] - integrals[k])/fabs(referenceIntegrals[k]) > 1.0e-6) &&
                                (errcount < 10))
                            {
                                printf ("2 ERROR: %"PRIu32" %"PRIu32" %"PRIu32" %"PRIu32": %le %le: %le\n",
                                    shellIndexM, shellIndexN, shellIndexP, shellIndexQ, referenceIntegrals[k], integrals[k],
                                    fabs(referenceIntegrals[k] - integrals[k])/fabs(referenceIntegrals[k]));
                                errcount++;
                            }

                        }
                    }   
                } else {
                    printf ("ERROR: nints0 %d nints %d\n", referenceIntegralsCount, integralsCount);
                }
                if (errcount > 0) {
                    goto end;
                }
            }
            
            free(shellIndicesN);
            free(shellIndicesQ);

        }
    }

    printf("Done\n\n");
    const double totalTime = (double)totalTicks / (double)getTimerFrequency();
    printf("Number of calls: %"PRIu64", Number of integrals: %"PRIu64"\n", totalCalls, totalIntegralsCount);
    printf("Total time: %.4lf secs\n", totalTime);
    printf("Average time per call: %.3lf us\n", totalTime / totalCalls * 1.0e+6);

end:
    free(referenceIntegrals);
    CInt_destroyERD(erd);
    CInt_destroyBasisSet(basis);
    free(shellptr);
    free(shellid);
    free(shellvalue);
    free(shellrid);
    
    fclose(ref_data_file);
    return 0;
}
