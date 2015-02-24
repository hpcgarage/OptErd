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
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <omp.h>

#include "CInt.h"
#include "screening.h"

#define CHUNK_SIZE 64
#define TASK_SIZE 5
#define CPU_ONLY
//#define MIC_ONLY

#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define FREE  alloc_if(0) free_if(1)

#pragma offload_attribute(push, target(mic))
ERD_t erd_mic;
double *D;

static void update_F (BasisSet_t basis_,
        int M, int N, int P, int Q,
        double *integrals,
        double *F, int ldm)
{
    int startM;
    int startN;
    int startP;
    int startQ;
    int dimM;
    int dimN;
    int dimP;
    int dimQ;
    int iM;
    int iN;
    int iP;
    int iQ;
    double ivalue;
    int fM;
    int fN;
    int fP;
    int fQ;

    startM = CInt_getFuncStartInd (basis_, M);
    startN = CInt_getFuncStartInd (basis_, N);
    startP = CInt_getFuncStartInd (basis_, P);
    startQ = CInt_getFuncStartInd (basis_, Q);
    dimM = CInt_getShellDim (basis_, M);
    dimN = CInt_getShellDim (basis_, N);
    dimP = CInt_getShellDim (basis_, P);
    dimQ = CInt_getShellDim (basis_, Q);		   
    for (iM = 0; iM < dimM; iM++)
    {
        for (iN = 0; iN < dimN; iN++)
        {
            for (iP = 0; iP < dimP; iP++)
            {
                for (iQ = 0; iQ < dimQ; iQ++)
                {
                    ivalue = integrals[iM + dimM * (iN + dimN * (iP + dimP * iQ))];
                    fM = startM + iM;
                    fN = startN + iN;
                    fP = startP + iP;
                    fQ = startQ + iQ;                    
                    F[fM * ldm + fN] += 0.5 * ivalue * D[fP * ldm + fQ];
                    F[fP * ldm + fQ] += 0.5 * ivalue * D[fM * ldm + fN];
                    F[fM * ldm + fP] += ivalue * D[fN * ldm + fQ];                                        
                }
            }
        }
    }
}



void processTasks(BasisSet_t basis_, ERD_t erd_, double *F,
        int *taskArray, int numTasks, int ldm, int sub_mat_size,
        int sub_mat_size_aligned, int nthreads)
{
    int i;

#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        int nints;
        int totalInts = 0;
        double *my_F = F + tid * sub_mat_size_aligned;

#pragma omp for schedule(dynamic, CHUNK_SIZE)
        for(i = 0; i < numTasks; i++)
        {
            int M = taskArray[i * TASK_SIZE + 0];
            int N = taskArray[i * TASK_SIZE + 1];
            int P = taskArray[i * TASK_SIZE + 2];
            int Q = taskArray[i * TASK_SIZE + 3];
            double *ints;
            CInt_computeShellQuartet (basis_, erd_, tid, M, N, P, Q, &ints, &nints);
            totalInts += nints;
            if (nints > 0)
                update_F (basis_, M, N, P, Q, ints, my_F, ldm);
        }
#pragma omp for
        for(i = 0; i < sub_mat_size; i++)
        {
            int id;
            for(id = 1; id < nthreads; id++)
            {
                F[i] += F[id * sub_mat_size_aligned + i];
            }
        }
    }
}

#pragma offload_attribute(pop)

void runTasksInHeteroMode(BasisSet_t basis, ERD_t erd, int *taskArray,
        int numTasks, double *F_hetero, double *F_hetero_cpu,
        int numTasks_mic, int ldm,
        int sub_mat_size, int sub_mat_size_aligned, int nthreads, int nthreads_mic,
        int num_devices)
{
    int i;
    int mic_id;
    int tasks_done_so_far = 0;
    int exp_ints_so_far = 0;
    long clock_start, clock_end;
    printf("runTasksInHeteroMode\n");


    // Initialize MIC cards for current set of tasks
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
        clock_start = __rdtsc();
#pragma offload target(mic:mic_id) \
        in(nthreads_mic, sub_mat_size_aligned) \
        in(D: length(sub_mat_size) REUSE) \
        in(F_hetero: length(0) REUSE)
        {
            int tid;
#pragma omp parallel for num_threads(nthreads_mic)
            for(tid = 0; tid < nthreads_mic; tid++)
            {
                memset(F_hetero + i * sub_mat_size_aligned, 0, sizeof(double) * sub_mat_size_aligned);
            }
        }
        clock_end = __rdtsc();
        printf("Ticks on hetero initializing mic:%d = %lf\n", mic_id, ((double)(clock_end - clock_start)) / 1000000000); 
    }
    
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
        int expectedTotalInts_mic = taskArray[(tasks_done_so_far + numTasks_mic) * 5 + 4] - exp_ints_so_far;

        printf("mic_id = %d, numTasks_mic = %d, expectedTotalInts_mic=%d, tasks_done_so_far = %d\n", mic_id, numTasks_mic, expectedTotalInts_mic, tasks_done_so_far);

#pragma offload target(mic:mic_id) \
        in(nthreads_mic, sub_mat_size, sub_mat_size_aligned, ldm) \
        in(numTasks_mic) \
        nocopy(D: length(sub_mat_size) alloc_if(0) free_if(0)) \
        in(taskArray[tasks_done_so_far * 5: numTasks_mic * 5]: into(taskArray[0: numTasks_mic * 5]) alloc_if(0) free_if(0)) \
        nocopy(basis_mic, erd_mic) \
        out(F_hetero[0 : sub_mat_size]: into(F_hetero[mic_id * sub_mat_size_aligned : sub_mat_size]) alloc_if(0) free_if(0)) \
        signal(F_hetero)
        processTasks(basis_mic, erd_mic, F_hetero,
                taskArray, numTasks_mic, ldm, sub_mat_size, 
                sub_mat_size_aligned, nthreads_mic);

        tasks_done_so_far += numTasks_mic;
        exp_ints_so_far = taskArray[numTasks_mic * 5 + 4];
    }
    printf("launched tasks on MIC. tasks_done_so_far = %d\n", tasks_done_so_far);

    // Process remaining tasks on CPU 
    int *taskArray_cpu = taskArray + tasks_done_so_far * TASK_SIZE;
    int numTasks_cpu = numTasks - tasks_done_so_far;
    processTasks(basis, erd, F_hetero_cpu,
            taskArray_cpu, numTasks_cpu, ldm, sub_mat_size,
            sub_mat_size_aligned, nthreads);
    printf("processed tasks on CPU\n");
#if 1
    if(numTasks_mic > 0)
    {
        for(mic_id = 0; mic_id < num_devices; mic_id++)
        {
#pragma offload target(mic:mic_id) wait(F_hetero)
            {
            }
        }
    }
#endif
    printf("tasks on MIC finished\n");
    for(i = 0; i < sub_mat_size; i++)
    {
        int id;
        for(id = 1; id < num_devices; id++)
        {
            F_hetero[i] 
                += F_hetero[id * sub_mat_size_aligned + i];
        }
        F_hetero[i] += F_hetero_cpu[i];
    }
}



int main (int argc, char **argv)
{
    int ns;
    int nnz;
    double totalcalls = 0;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    //struct timeval tv1;
    //struct timeval tv2;
    //double timepass;   
    double totalintls = 0;
    BasisSet_t basis;
    ERD_t erd;
    //double *D;
    double *F_hetero;
    int errcount;
    int ldm;
    int i;
    int j;
    int start1;
    int end1;
    int start2;
    int end2;
    double value1;
    double value2;
    //int k;
    double norm;
    int *taskArray;

    if (argc != 6)
    {
        printf ("Usage: %s <basisset> <xyz> <nthreads> <nthreads_mic> <mic_fraction>\n", argv[0]);
        return 0;
    }
    errcount = 0;
    int nthreads = atoi (argv[3]);
    omp_set_num_threads (nthreads);
    int nthreads_mic = atoi (argv[4]);
    // load basis set   
    CInt_createBasisSet (&basis);
    CInt_loadBasisSet (basis, argv[1], argv[2]);
    schwartz_screening (basis, &shellptr, &shellid, &shellrid, &shellvalue, &nnz);

    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));
    printf ("  nthreads\t= %d\n", nthreads);

    ldm = CInt_getNumFuncs (basis);
    int sub_mat_size = ldm * ldm;
    int sub_mat_size_aligned = ((sub_mat_size + 63)/64) * 64;
    D = (double *)malloc (sizeof(double) * sub_mat_size);
    assert (D != NULL);
    for (i = 0; i < ldm * ldm; i++)
    {
        D[i] = 1.0;
    }

    // Find the number of MIC cards available
    int num_devices = 0;   
#ifdef __INTEL_OFFLOAD
    num_devices = _Offload_number_of_devices();
#endif

    if(num_devices == 0)
    {
        printf("No target devices available. Exiting\n");
        exit(0);
    }
    else
    {
        printf("Number of Target devices installed: %d\n\n",num_devices);
    }

    printf ("Computing integrals ...\n");
    ns = CInt_getNumShells (basis);
    int mic_id;

    // Create ERD on MIC cards
    for (mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload target(mic:mic_id) nocopy(basis_mic, erd_mic) in(nthreads_mic)
        {
            CInt_createERD (basis_mic, &erd_mic, nthreads_mic);
#pragma omp parallel num_threads(nthreads_mic)
            {
                int thread_count = omp_get_thread_num();
            }
            /*erd_mic = (ERD_t *) malloc (sizeof (ERD_t) * nthreads_mic);
            assert (erd_mic != NULL);

#pragma omp parallel for num_threads(nthreads_mic)
            for(int i = 0; i < nthreads_mic; i++)
            {
                CInt_createERD (basis_mic, &(erd_mic[i]));
            }*/
        }
    }

    // Create ERD on CPU
    CInt_createERD (basis, &erd, nthreads);
    /*
    erd = (ERD_t *) malloc (sizeof (ERD_t) * nthreads);
    assert (erd != NULL);
#pragma omp parallel for
    for (i = 0; i < nthreads; i++)
    {
        CInt_createERD (basis, &(erd[i]));
    }
*/
    // Prepare all the tasks
    taskArray = (int *)_mm_malloc(ns * ns * ns * ns * 5 * sizeof(int), 2 * 1024 * 1024);

    int numTasks = 0;
    int expectedTotalInts = 0;
    long clock_start, clock_end;
    clock_start = __rdtsc();
    for (int M = 0; M < ns; M++)
    {
        start1 = shellptr[M];
        end1 = shellptr[M + 1];       
        for (i = start1; i < end1; i++)
        {
            int N = shellid[i];
            value1 = shellvalue[i];
            for (int P = 0; P < ns; P++)
            {
                start2 = shellptr[P];
                end2 = shellptr[P + 1];              
                for (j = start2; j < end2; j++)
                {
                    int Q = shellid[j];
                    value2 = shellvalue[j];
                    if (M > N || P > Q || (M + N) > (P + Q))
                        continue;
                    if (fabs(value1 * value2) < TOLSRC * TOLSRC)
                        continue;        
                    taskArray[numTasks * 5 + 0] = M;                
                    taskArray[numTasks * 5 + 1] = N;                
                    taskArray[numTasks * 5 + 2] = P;                
                    taskArray[numTasks * 5 + 3] = Q; 
                    taskArray[numTasks * 5 + 4] = expectedTotalInts;
                    numTasks++;
                    int dimM = CInt_getShellDim (basis, M);
                    int dimN = CInt_getShellDim (basis, N);
                    int dimP = CInt_getShellDim (basis, P);
                    int dimQ = CInt_getShellDim (basis, Q);		   
                    expectedTotalInts += dimM * dimN * dimP * dimQ;
                } /* for (j = start2; j < end2; j++) */
            } /* for (P = 0; P < ns; P++) */
        } /* for (i = start1; i < end1; i++) */
    } /* for (M = 0; M < ns; M++) */
    clock_end = __rdtsc();
    printf("Ticks for initializing taskArray = %lf\n", ((double)(clock_end - clock_start)) / 1000000000); 
    printf("numTasks = %d, ns * ns * ns * ns = %d\n", numTasks, ns * ns * ns * ns); 
    printf("expectedTotalInts = %d\n", expectedTotalInts);

#ifdef MIC_ONLY
// Run the tasks on MIC only
    double *F_mic = (double *)_mm_malloc (sizeof(double) * nthreads_mic * sub_mat_size_aligned, 2 * 1024 * 1024);
    assert (F_mic != NULL);
    memset (F_mic, 0, sizeof(double) * sub_mat_size_aligned);

    printf("Running on MIC: numTasks = %d, expectedTotalInts=%d\n", numTasks, expectedTotalInts);

// Allocate memory and send input data
#pragma offload target(mic:0) \
    in(nthreads_mic, sub_mat_size_aligned, ldm) \
    in(D: length(sub_mat_size) ALLOC) \
    in(taskArray: length(numTasks * TASK_SIZE) ALLOC) \
    nocopy(F_mic: length(nthreads_mic * sub_mat_size_aligned) ALLOC)
    {   
        memset(F_mic, 0, sizeof(double) * nthreads_mic * sub_mat_size_aligned);
    }

// Launch kernel on MIC
    clock_start = __rdtsc();
#pragma offload target(mic:0) \
    in(numTasks, nthreads_mic, sub_mat_size, ldm) \
    nocopy(D: REUSE) \
    nocopy(taskArray: REUSE) \
    nocopy(basis_mic, erd_mic: REUSE) \
    out(F_mic: length(sub_mat_size) REUSE)
    processTasks(basis_mic, erd_mic, F_mic,
            taskArray, numTasks, ldm, sub_mat_size,
            sub_mat_size_aligned, nthreads_mic);

    clock_end = __rdtsc();
    printf("Ticks on MIC = %lf\n", ((double)(clock_end - clock_start)) / 1000000000); 

// Free allocated buffers
#pragma offload_transfer target(mic:0) \
    nocopy(taskArray: length(0) FREE) \
    nocopy(D: length(0) FREE) \
    nocopy(F_mic: length(0) FREE)

#endif
#ifdef CPU_ONLY
    // Run the tasks on CPU only
    double *F = (double *)_mm_malloc (sizeof(double) * sub_mat_size_aligned * nthreads, 64);
    assert (F != NULL);
    memset (F, 0, sizeof(double) * sub_mat_size_aligned * nthreads);
    clock_start = __rdtsc();
    processTasks(basis, erd, F,
            taskArray, numTasks, ldm, sub_mat_size,
            sub_mat_size_aligned, nthreads);
    clock_end = __rdtsc();
    printf("Ticks on CPU = %lf\n", ((double)(clock_end - clock_start)) / 1000000000); 
#endif

    // hetero code
    // Allocate required arrays
    
    printf("sub_mat_size = %d, sub_mat_size_aligned = %d\n", sub_mat_size, sub_mat_size_aligned);

    F_hetero = (double *)_mm_malloc (sizeof(double) * sub_mat_size_aligned * (num_devices + 1), 64);
    assert (F_hetero != NULL);
    memset (F_hetero, 0, sizeof(double) * sub_mat_size_aligned * (num_devices + 1));

    double *F_hetero_cpu = (double *)_mm_malloc(sizeof(double) * sub_mat_size_aligned * nthreads, 64);
    assert (F_hetero_cpu != NULL);
    memset(F_hetero_cpu, 0, sizeof(double) * sub_mat_size_aligned * nthreads);

    double mic_fraction = atof(argv[5]);

    int numTasks_mic = numTasks * mic_fraction;
    // Ensure that the size of taskArray to be sent to mic is a multiple of 64 bytes.
    numTasks_mic = (numTasks_mic /16 ) * 16;

    printf("mic_fraction = %lf, numTasks_mic = %d\n", mic_fraction, numTasks_mic);
    // Allocate arrays on MIC cards
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
        clock_start = __rdtsc();
#pragma offload target(mic:mic_id) \
        in(nthreads_mic, sub_mat_size_aligned) \
        in(D: length(sub_mat_size) ALLOC) \
        nocopy(taskArray: length(numTasks_mic * 5) ALLOC) \
        nocopy(F_hetero: length(nthreads_mic * sub_mat_size_aligned) ALLOC)
        {
        }
        clock_end = __rdtsc();
        printf("Ticks on hetero allocating arrays on mic:%d = %lf\n", mic_id, ((double)(clock_end - clock_start)) / 1000000000); 
    }
    printf("allocated memory\n");


    long hetero_start, hetero_end;
    hetero_start = __rdtsc();
    runTasksInHeteroMode(basis, erd, taskArray, numTasks,
            F_hetero, F_hetero_cpu, numTasks_mic, ldm,
            sub_mat_size, sub_mat_size_aligned, nthreads, nthreads_mic, num_devices);
    hetero_end = __rdtsc();
    printf("Ticks on hetero = %lf\n", ((double)(hetero_end - hetero_start)) / 1000000000); 

    // Free arrays on MIC
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload target(mic:mic_id) \
        nocopy(D: length(0) FREE) \
        nocopy(F_hetero: length(0) FREE)
        {
        }
    }

    //free(taskArray);

    printf ("Done\n");
    printf ("\n");
    printf ("Number of calls: %.6le, Number of integrals: %.6le\n",
            totalcalls, totalintls);
    //printf ("Total time: %.4lf secs\n", timepass);
    //printf ("Average time per call: %.3le us\n",
    //        1000.0 * 1000.0 * timepass / totalcalls);

#ifdef CPU_ONLY
    cblas_daxpy (ldm * ldm, -1.0, F_hetero, 1, F, 1); 
    norm = dlange ("F", &ldm, &ldm, F, &ldm, NULL);
    printf ("F diff is %le\n", norm);
    if (norm > 1e-6)
    {
        printf ("ERROR: F diff is larger than 1e-6\n");
    }
#endif

    for (mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload target(mic:mic_id) nocopy(erd_mic)
        {
            CInt_destroyERD (erd_mic);
            /*
            for(int i = 0; i < nthreads_mic; i++)
            {
                CInt_destroyERD (erd_mic[i]);
            }
            free(erd_mic);
            */
        }
    }
    CInt_destroyERD (erd);
    /*
    for (i = 0; i < nthreads; i++)
    {
        CInt_destroyERD (erd[i]);
    }
    free (erd);
    */
    CInt_destroyBasisSet (basis);

    free (shellptr);
    free (shellid);
    free (shellvalue);
    free (shellrid);
    free (D);
#ifdef CPU_ONLY
    _mm_free (F);
#endif
#ifdef MIC_ONLY
    _mm_free (F_mic);
#endif
    _mm_free (F_hetero);

    return 0;
}
