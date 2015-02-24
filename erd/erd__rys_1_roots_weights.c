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

#include "erd.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__RYS_1_ROOTS_WEIGHTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation returns Rys polynomial roots and weights */
/*                in case the number of roots and weights required */
/*                is = 1. All T's are treated at once so the complete */
/*                set of roots and weights is returned. */
/*                For the moment taken essentially unchanged from the */
/*                GAMESS package (routine RTS123, but removing their */
/*                'spaghetti' code from the 70's of unreadable */
/*                internested IFs and GOTOs!). */
/*                One interesting aspect of the GAMESS routines is that */
/*                their code returns scaled roots, i.e. their roots */
/*                do not ly between the range 0 and 1. To get to the */
/*                proper roots as needed for our package, we simply */
/*                set: */
/*                   root (our) = root (gamess) / (1 + root (games)) */
/*                  Input: */
/*                    NT           =  # of T-exponents */
/*                    TVAL         =  the set of NT T-exponents defining */
/*                                    the Rys weight functions */
/*                  Output: */
/*                    RTS          =  all NT quadrature roots */
/*                    WTS          =  all NT quadrature weights */
/* ------------------------------------------------------------------------ */
void erd__rys_1_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]) {
    int jump1[34] =
        { 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7
    };

    double e;
    int n;
    double t, x, f1, r1, w1;
    int tcase;

/* ------------------------------------------------------------------------ */
/*                 ******************************** */
/*             ... *  # of roots and weights = 1  * */
/*                 ******************************** */
    for (n = 0; n < nt; ++n)
    {
        t = tval[n];
        if (t <= 3e-7)
        {
/*             ...T-range: T essentially 0 */
            r1 = .5 - t * .2;
            wts[n] *= 1. - t * .333333333333333;
            rts[n] = r1 / (r1 + 1.);
            goto L100;
        }
        tcase = (int) MIN ((t + 1.0), 34.);
        switch (jump1[tcase - 1])
        {
        case 1:
            goto L1100;
        case 2:
            goto L1200;
        case 3:
            goto L1300;
        case 4:
            goto L1400;
        case 5:
            goto L1500;
        case 6:
            goto L1600;
        case 7:
            goto L1700;
        }

/*             ...T-range: 0 < T < 1 */
      L1100:
        f1 = ((((((((t * -8.36313918003957e-8 + 1.21222603512827e-6) * t -
                    1.15662609053481e-5) * t + 9.25197374512647e-5) * t -
                  6.40994113129432e-4) * t + .00378787044215009) * t -
                .0185185172458485) * t + .0714285713298222) * t -
              .199999999997023) * t + .333333333333318;
        w1 = (t + t) * f1 + exp (-t);
        r1 = f1 / (w1 - f1);
        wts[n] *= w1;
        rts[n] = r1 / (r1 + 1.);
        goto L100;


/*             ...T-range: 1 =< T < 3 */
      L1200:
        x = t - 2.;
        f1 = ((((((((((x * -1.61702782425558e-10 + 1.96215250865776e-9) * x -
                      2.14234468198419e-8) * x + 2.17216556336318e-7) * x -
                    1.98850171329371e-6) * x + 1.62429321438911e-5) * x -
                  1.16740298039895e-4) * x + 7.24888732052332e-4) * x -
                .00379490003707156) * x + .0161723488664661) * x -
              .0529428148329736) * x + .115702180856167;
        w1 = (t + t) * f1 + exp (-t);
        r1 = f1 / (w1 - f1);
        wts[n] *= w1;
        rts[n] = r1 / (r1 + 1.);
        goto L100;


/*             ...T-range: 3 =< T < 5 */
      L1300:
        x = t - 4.;
        f1 = ((((((((((x * -2.62453564772299e-11 + 3.24031041623823e-10) * x
                      - 3.614965656163e-9) * x + 3.760256799971e-8) * x -
                    3.553558319675e-7) * x + 3.022556449731e-6) * x -
                  2.290098979647e-5) * x + 1.526537461148e-4) * x -
                8.81947375894379e-4) * x + .00433207949514611) * x -
              .0175257821619926) * x + .0528406320615584;
        w1 = (t + t) * f1 + exp (-t);
        r1 = f1 / (w1 - f1);
        wts[n] *= w1;
        rts[n] = r1 / (r1 + 1.);
        goto L100;


/*             ...T-range: 5 =< T < 10 */
      L1400:
        e = exp (-t);
        x = 1. / t;
        w1 = ((((((x * .46897511375022 - .69955602298985) * x +
                  .53689283271887) * x - .32883030418398) * x +
                .24645596956002) * x - .49984072848436) * x -
              3.1501078774085e-6) * e + sqrt (x * .785398163397448);
        f1 = (w1 - e) / (t + t);
        r1 = f1 / (w1 - f1);
        wts[n] *= w1;
        rts[n] = r1 / (r1 + 1.);
        goto L100;


/*             ...T-range: 10 =< T < 15 */
      L1500:
        e = exp (-t);
        x = 1. / t;
        w1 = (((x * -.18784686463512 + .22991849164985) * x - .49893752514047)
              * x - 2.1916512131607e-5) * e + sqrt (x * .785398163397448);
        f1 = (w1 - e) / (t + t);
        r1 = f1 / (w1 - f1);
        wts[n] *= w1;
        rts[n] = r1 / (r1 + 1.);
        goto L100;


/*             ...T-range: 15 =< T < 33 */
      L1600:
        e = exp (-t);
        x = 1. / t;
        w1 = ((x * .1962326414943 - .4969524146449) * x - 6.0156581186481e-5)
            * e + sqrt (x * .785398163397448);
        f1 = (w1 - e) / (t + t);
        r1 = f1 / (w1 - f1);
        wts[n] *= w1;
        rts[n] = r1 / (r1 + 1.);
        goto L100;


/*             ...T-range: T >= 33 */
      L1700:
        wts[n] *= sqrt (.785398163397448 / t);
        rts[n] = .5 / t;
      L100:
        ;
    }
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
