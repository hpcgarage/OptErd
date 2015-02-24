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

#include "CInt.h"
#include "screening.h"


#pragma offload_attribute(push, target(mic))
ERD_t erd_mic;
double integrals_mic[10000];
#pragma offload_attribute(pop)


static void update_F (BasisSet_t basis,
                      int M, int N, int P, int Q,
                      double *integrals,
                      double *D, double *F, int ldm)
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
    
    startM = CInt_getFuncStartInd (basis, M);
    startN = CInt_getFuncStartInd (basis, N);
    startP = CInt_getFuncStartInd (basis, P);
    startQ = CInt_getFuncStartInd (basis, Q);
    dimM = CInt_getShellDim (basis, M);
    dimN = CInt_getShellDim (basis, N);
    dimP = CInt_getShellDim (basis, P);
    dimQ = CInt_getShellDim (basis, Q);		   
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


int main (int argc, char **argv)
{
    int M;
    int N;
    int P;
    int Q;
    int nints_mic;
    int ns;
    int nnz;
    double totalcalls = 0;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    struct timeval tv1;
    struct timeval tv2;
    double timepass;   
    double totalintls = 0;
    BasisSet_t basis;
    ERD_t erd;
    double *D;
    double *F_mic;
    double *F;
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
    int nints;
    int k;
    double *integrals;
    double norm;
    
    if (argc != 3)
    {
        printf ("Usage: %s <basisset> <xyz>\n", argv[0]);
        return 0;
    }
    errcount = 0;
    // load basis set   
    CInt_createBasisSet (&basis);
    CInt_loadBasisSet (basis, argv[1], argv[2]);
    schwartz_screening (basis, &shellptr, &shellid, &shellrid, &shellvalue, &nnz);

    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));

    ldm = CInt_getNumFuncs (basis);
    D = (double *)malloc (sizeof(double) * ldm * ldm);
    F_mic = (double *)malloc (sizeof(double) * ldm * ldm);
    F = (double *)malloc (sizeof(double) * ldm * ldm);
    assert (D != NULL &&
            F != NULL &&
            F_mic != NULL);
    memset (F_mic, 0, sizeof(double) * ldm *ldm);
    memset (F, 0, sizeof(double) * ldm *ldm);
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
    timepass = 0.0;
    for (i = 0; i < num_devices; i++)
    {
        #pragma offload target(mic:i) nocopy(basis_mic, erd_mic)
        {
            CInt_createERD (basis_mic, &erd_mic, 1);
        }
    }
    CInt_createERD (basis, &erd, 1);
    
    for (M = 0; M < ns; M++)
    {
        start1 = shellptr[M];
        end1 = shellptr[M + 1];       
        for (i = start1; i < end1; i++)
        {
            N = shellid[i];
            value1 = shellvalue[i];
            for (P = 0; P < ns; P++)
            {
                start2 = shellptr[P];
                end2 = shellptr[P + 1];              
                for (j = start2; j < end2; j++)
                {
                    Q = shellid[j];
                    value2 = shellvalue[j];
                    if (M > N || P > Q || (M + N) > (P + Q))
                        continue;
                    if (fabs(value1 * value2) < TOLSRC * TOLSRC)
                        continue;                        
                    totalcalls = totalcalls + 1;
                    gettimeofday (&tv1, NULL);
                    #pragma offload target(mic:0) \
                            in(M, N, P, Q) \
                            nocopy(basis_mic, erd_mic) \
                            out(integrals_mic, nints_mic)
                    {
                        double *ints;
                        CInt_computeShellQuartet (basis_mic, erd_mic, 0, M, N, P, Q,
                                &ints, &nints_mic);
                        memcpy(integrals_mic, ints, nints_mic * sizeof(double));
                    }
                    gettimeofday (&tv2, NULL);
                    timepass += (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
                    totalintls = totalintls + nints_mic;
                    // CPU call
                    CInt_computeShellQuartet (basis, erd, 0, M, N, P, Q, &integrals, &nints);
                    
                    // compare results
                    if (nints_mic == 0 && nints == 0)
                    {
                        continue;
                    }
                    else if (nints_mic != 0 && nints == 0)
                    {
                        for (k = 0; k < nints_mic; k++)
                        {
                            if (integrals_mic[k] > 1e-10)
                            {
                                printf ("ERROR: %d %d %d %d: %le %le\n",
                                    M, N, P, Q, 0.0, integrals_mic[k]);
                                errcount++;
                            }
                        }
                    }
                    else if (nints_mic == 0 && nints != 0)
                    {
                        for (k = 0; k < nints; k++)
                        {
                            if (integrals[k] > 1e-10)
                            {
                                printf ("ERROR: %d %d %d %d: %le %le\n",
                                    M, N, P, Q, integrals[k], 0.0);
                                errcount++;
                            }
                        }
                    
                    }
                    else if (nints == nints_mic && nints_mic != 0)
                    {
                     
                        for (k = 0; k < nints_mic; k++)
                        {
                            if (fabs(integrals[k]) < 1e-6 ||
                                fabs(integrals_mic[k]) < 1e-6)
                            {
                                if (fabs(integrals[k] - integrals_mic[k]) > 1e-10)
                                {
                                    printf ("1 ERROR: %d %d %d %d: %le %le\n",
                                        M, N, P, Q, integrals[k], integrals_mic[k]);
                                    errcount++;
                                }
                            }
                            else
                            {
                                if (fabs(integrals[k] - integrals_mic[k])/fabs(integrals[k]) >
                                    1e-6 && errcount < 10)
                                {
                                    printf ("2 ERROR: %d %d %d %d: %le %le: %le\n",
                                        M, N, P, Q, integrals[k], integrals_mic[k],
                                        fabs(integrals[k] - integrals_mic[k])/fabs(integrals[k]));
                                    errcount++;
                                }

                            }
                        }   
                    }
                    else
                    {
                        printf ("ERROR: nints0 %d nints %d\n", nints_mic, nints);
                    }

                    if (errcount > 10)
                    {
                        goto end;
                    }
                    update_F (basis, M, N, P, Q, integrals, D, F, ldm);
                    update_F (basis, M, N, P, Q, integrals_mic, D, F_mic, ldm);
                } /* for (j = start2; j < end2; j++) */
            } /* for (P = 0; P < ns; P++) */
        } /* for (i = start1; i < end1; i++) */
    } /* for (M = 0; M < ns; M++) */
    
    printf ("Done\n");
    printf ("\n");
    printf ("Number of calls: %.6le, Number of integrals: %.6le\n",
            totalcalls, totalintls);
    printf ("Total time: %.4lf secs\n", timepass);
    printf ("Average time per call: %.3le us\n",
            1000.0 * 1000.0 * timepass / totalcalls);

    cblas_daxpy (ldm * ldm, -1.0, F_mic, 1, F, 1); 
    norm = dlange ("F", &ldm, &ldm, F, &ldm, NULL);
    printf ("F diff is %le\n", norm);
    if (norm > 1e-6)
    {
        printf ("ERROR: F diff is larger than 1e-6\n");
    }
end:
    for (i = 0; i < num_devices; i++)
    {
        #pragma offload target(mic:i) nocopy(erd_mic)
        {
            CInt_destroyERD (erd_mic);
        }
    }
    CInt_destroyERD (erd);
    CInt_destroyBasisSet (basis);
    free (shellptr);
    free (shellid);
    free (shellvalue);
    free (shellrid);
    free (D);
    free (F);
    free (F_mic);
    
    return 0;
}
