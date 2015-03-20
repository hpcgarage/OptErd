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
#include "erdutil.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__RYS_4_ROOTS_WEIGHTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation returns Rys polynomial roots and weights */
/*                in case the number of roots and weights required */
/*                is = 4. All T's are treated at once so the complete */
/*                set of roots and weights is returned. */
/*                For the moment taken essentially unchanged from the */
/*                GAMESS package (routine ROOT4, but removing their */
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
/*                    NTGQP        =  # of roots times # of T-exponents */
/*                                    (= 4 * NT) */
/*                    TVAL         =  the set of NT T-exponents defining */
/*                                    the Rys weight functions */
/*                  Output: */
/*                    RTS          =  all NTGQP quadrature roots */
/*                    WTS          =  all NTGQP quadrature weights */
/* ------------------------------------------------------------------------ */
void erd__rys_4_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]) {
    int jump4[54] =
        { 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 8
    };

    double e;
    int m, n;
    double t, x, y, r1, r2, r3, r4, w1, w2, w3, w4;
    int tcase;

/* ------------------------------------------------------------------------ */
/*                 ******************************** */
/*             ... *  # of roots and weights = 4  * */
/*                 ******************************** */
    m = 0;
    for (n = 0; n < nt; ++n)
    {
        t = tval[n];
        if (t <= 3e-7)
        {
/*             ...T-range: T essentially 0 */
            r1 = .0348198973061471 - t * .00409645850660395;
            r2 = .381567185080042 - t * .0448902570656719;
            r3 = 1.73730726945891 - t * .204389090547327;
            r4 = 11.8463056481549 - t * 1.39368301742312;
            wts[m] *= .362683783378362 - t * .0313844305713928;
            wts[m + 1] *= .313706645877886 - t * .0898046242557724;
            wts[m + 2] *= .222381034453372 - t * .129314370958973;
            wts[m + 3] *= .101228536290376 - t * .0828299075414321;
            rts[m] = r1 / (r1 + 1.);
            rts[m + 1] = r2 / (r2 + 1.);
            rts[m + 2] = r3 / (r3 + 1.);
            rts[m + 3] = r4 / (r4 + 1.);
            m += 4;
            goto L400;
        }

        tcase = (int) MIN ((t + 1.0), 54.);
        switch (jump4[tcase - 1])
        {
        case 1:
            goto L4100;
        case 2:
            goto L4200;
        case 3:
            goto L4300;
        case 4:
            goto L4400;
        case 5:
            goto L4500;
        case 6:
            goto L4600;
        case 7:
            goto L4700;
        case 8:
            goto L4800;
        }


/*             ...T-range: 0 < T < 1 */
      L4100:
        wts[m] *= ((((((t * -1.14649303201279e-8 + 1.88015570196787e-7) * t -
                      2.33305875372323e-6) * t + 2.68880044371597e-5) * t -
                    2.94268428977387e-4) * t + .00306548909776613) * t -
                  .0313844305680096) * t + .362683783378335;
        wts[m + 1] *= ((((((((t * -4.11720483772634e-9 + 6.54963481852134e-8) *
                            t - 7.20045285129626e-7) * t +
                           6.93779646721723e-6) * t -
                          6.05367572016373e-5) * t +
                         4.74241566251899e-4) * t - .00326956188125316) * t +
                       .0191883866626681) * t - .0898046242565811) * t +
            .313706645877886;
        wts[m + 2] *=
            ((((((((t * -3.41688436990215e-8 + 5.07238960340773e-7) * t -
                   5.0167562840822e-6) * t + 4.20363420922845e-5) * t -
                 3.08040221166823e-4) * t + .00194431864731239) * t -
               .0102477820460278) * t + .0428670143840073) * t -
             .129314370962569) * t + .222381034453369;
        wts[m + 3] *=
            (((((((((t * 4.99660550769508e-9 - 7.9458596331012e-8) * t +
                    8.359072409485e-7) * t - 7.42236921061e-6) * t +
                  5.76337430816e-5) * t - 3.86645606718233e-4) * t +
                .00218417516259781) * t - .00999791027771119) * t +
              .034879109737737) * t - .0828299075413889) * t +
            .101228536290376;
        r1 = ((((((t * -1.95309614628539e-10 + 5.19765728707592e-9) * t -
                  1.01756452250573e-7) * t + 1.72365935872131e-6) * t -
                2.61203523522184e-5) * t + 3.5292130876988e-4) * t -
              .00409645850658433) * t + .0348198973061469;
        r2 = (((((t * -1.89554881382342e-8 + 3.07583114342365e-7) * t +
                 1.270981734393e-6) * t - 1.417298563884e-4) * t +
               .003226979163176) * t - .0448902570678178) * t +
            .381567185080039;
        r3 = ((((((t * 1.77280535300416e-9 + 3.36524958870615e-8) * t -
                  2.58341529013893e-7) * t - 1.1364489566232e-5) * t -
                7.91549618884063e-5) * t + .0103825827346828) * t -
              .204389090525137) * t + 1.73730726945889;
        r4 = (((((t * -5.61188882415248e-8 - 2.4948073307246e-7) * t +
                 3.428685057114e-6) * t + 1.679007454539e-4) * t +
               .04722855585715) * t - 1.39368301737828) * t +
            11.8463056481543;
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: 1 =< T < 5 */
      L4200:
        x = t - 3.;
        wts[m] *= ((((((((((x * -4.65801912689961e-14 + 7.586695071068e-13) *
                          x - 1.186387548048e-11) * x +
                         1.862334710665e-10) * x - 2.799399389539e-9) * x +
                       4.148972684255e-8) * x - 5.9335680796e-7) * x +
                     8.168349266115e-6) * x - 1.08989176177409e-4) * x +
                   .00141357961729531) * x - .0187588361833659) * x +
            .289898651436026;
        wts[m + 1] *=
            ((((((((((((x * -1.46345073267549e-14 +
                        2.25644205432182e-13) * x - 3.116258693847e-12) * x +
                      4.32190875661e-11) * x - 5.673270062669e-10) * x +
                    7.00629596296e-9) * x - 8.120186517e-8) * x +
                  8.77529464577e-7) * x - 8.77829235749024e-6) * x +
                8.04372147732379e-5) * x - 6.64149238804153e-4) * x +
              .00481181506827225) * x - .0288982669486183) * x +
            .156247249979288;
        wts[m + 2] *=
            (((((((((((((x * 9.06812118895365e-15 -
                         1.40541322766087e-13) * x + 1.919270015269e-12) * x -
                       2.60513573901e-11) * x + 3.299685839012e-10) * x -
                     3.86354139348735e-9) * x + 4.16265847927498e-8) * x -
                   4.0946283547147e-7) * x + 3.64018881086111e-6) * x -
                 2.88665153269386e-5) * x + 2.00515819789028e-4) * x -
               .00118791896897934) * x + .00575223633388589) * x -
             .0209400418772687) * x + .0485368861938873;
        wts[m + 3] *=
            ((((((((((((((x * -9.74835552342257e-16 +
                          1.57857099317175e-14) * x -
                         2.249993780112e-13) * x + 3.173422008953e-12) * x -
                       4.16115945968e-11) * x + 5.021343560166e-10) * x -
                     5.545047534808e-9) * x + 5.554146993491e-8) * x -
                   4.99048696190133e-7) * x + 3.96650392371311e-6) * x -
                 2.73816413291214e-5) * x + 1.60106988333186e-4) * x -
               7.64560567879592e-4) * x + .00281330044426892) * x -
             .00716227030134947) * x + .00966077262223353;
        r1 = (((((((((x * -1.48570633747284e-15 - 1.33273068108777e-13) * x +
                     4.06854369667e-12) * x - 9.163164161821e-11) * x +
                   2.046819017845e-9) * x - 4.03076426299031e-8) * x +
                 7.29407420660149e-7) * x - 1.23118059980833e-5) * x +
               1.88796581246938e-4) * x - .00253262912046853) * x +
            .0251198234505021;
        r2 = (((((((((x * 1.35830583483312e-13 - 2.29772605964836e-12) * x -
                     3.821500128045e-12) * x + 6.844424214735e-10) * x -
                   1.048063352259e-8) * x + 1.50083186233363e-8) * x +
                 3.48848942324454e-6) * x - 1.08694174399193e-4) * x +
               .00208048885251999) * x - .0291205805373793) * x +
            .272276489515713;
        r3 = (((((((((x * 5.02799392850289e-13 + 1.07461812944084e-11) * x -
                     1.482277886411e-10) * x - 2.153585661215e-9) * x +
                   3.654087802817e-8) * x + 5.1592957583012e-7) * x -
                 9.52388379435709e-6) * x - 2.16552440036426e-4) * x +
               .0090355146956832) * x - .145505469175613) * x +
            1.21449092319186;
        r4 = (((((((((x * -1.08510370291979e-12 + 6.41492397277798e-11) * x +
                     7.542387436125e-10) * x - 2.213111836647e-9) * x -
                   1.448228963549e-7) * x - 1.95670833237101e-6) * x -
                 1.07481314670844e-5) * x + 1.49335941252765e-4) * x +
               .0487791531990593) * x - 1.10559909038653) * x +
            8.0950202861178;
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: 5 =< T < 10 */
      L4300:
        x = t - 7.5;
        wts[m] *= ((((((((((x * -1.65995045235997e-15 + 6.91838935879598e-14) *
                          x - 9.131223418888e-13) * x +
                         1.403341829454e-11) * x - 3.672235069444e-10) * x +
                       6.36696254699e-9) * x - 1.039220021671e-7) * x +
                     1.959098751715e-6) * x - 3.33474893152939e-5) * x +
                   5.72164211151013e-4) * x - .0105583210553392) * x +
            .226696066029591;
        wts[m + 1] *=
            ((((((((((((x * -3.57248951192047e-16 +
                        6.25708409149331e-15) * x - 9.657033089714e-14) * x +
                      1.507864898748e-12) * x - 2.33252225611e-11) * x +
                    3.428545616603e-10) * x - 4.698730937661e-9) * x +
                  6.21997763513e-8) * x - 7.83008889613661e-7) * x +
                9.08621687041567e-6) * x - 9.86368311253873e-5) * x +
              9.69632496710088e-4) * x - .00814594214284187) * x +
            .0850218447733457;
        wts[m + 2] *=
            (((((((((((((x * 1.64742458534277e-16 - 2.6851226592841e-15) * x +
                        3.788890667676e-14) * x - 5.508918529823e-13) * x +
                      7.555896810069e-12) * x - 9.69039768312637e-11) * x +
                    1.16034263529672e-9) * x - 1.28771698573873e-8) * x +
                  1.31949431805798e-7) * x - 1.23673915616005e-6) * x +
                1.04189803544936e-5) * x - 7.79566003744742e-5) * x +
              5.03162624754434e-4) * x - .00255138844587555) * x +
            .0113250730954014;
        wts[m + 3] *=
            ((((((((((((((x * -1.55714130075679e-17 +
                          2.57193722698891e-16) * x -
                         3.626606654097e-15) * x + 5.234734676175e-14) * x -
                       7.067105402134e-13) * x + 8.79351266489e-12) * x -
                     1.006088923498e-10) * x + 1.050565098393e-9) * x -
                   9.91517881772662e-9) * x + 8.35835975882941e-8) * x -
                 6.19785782240693e-7) * x + 3.95841149373135e-6) * x -
               2.11366761402403e-5) * x + 9.00474771229507e-5) * x -
             2.78777909813289e-4) * x + 5.26543779837487e-4;
        r1 = (((((((((x * 4.64217329776215e-15 - 6.27892383644164e-15) * x +
                     3.462236347446e-13) * x - 2.92722935535e-11) * x +
                   5.090355371676e-10) * x - 9.97272656345253e-9) * x +
                 2.37835295639281e-7) * x - 4.60301761310921e-6) * x +
               8.42824204233222e-5) * x - .00137983082233081) * x +
            .0166630865869375;
        r2 = (((((((((x * 2.93981127919047e-14 + 8.47635639065744e-13) * x -
                     1.446314544774e-11) * x - 6.149155555753e-12) * x +
                   8.484275604612e-10) * x - 6.10898827887652e-8) * x +
                 2.39156093611106e-6) * x - 5.35837089462592e-5) * x +
               .00100967602595557) * x - .0157769317127372) * x +
            .174853819464285;
        r3 = ((((((((((x * 2.93523563363e-14 - 6.4004177666702e-14) * x -
                      2.695740446312e-12) * x + 1.027082960169e-10) * x -
                    5.82203865678e-10) * x - 3.159991002539e-8) * x +
                  4.327249251331e-7) * x + 4.856768455119e-6) * x -
                2.54617989427762e-4) * x + .00554843378106589) * x -
              .0795013029486684) * x + .720206142703162;
        r4 = (((((((((((x * -1.62212382394553e-14 +
                        7.68943641360593e-13) * x + 5.764015756615e-12) * x -
                      1.380635298784e-10) * x - 1.476849808675e-9) * x +
                    1.84347052385605e-8) * x + 3.34382940759405e-7) * x -
                  1.39428366421645e-6) * x - 7.50249313713996e-5) * x -
                6.26495899187507e-4) * x + .0469716410901162) * x -
              .666871297428209) * x + 4.11207530217806;
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: 10 =< T < 15 */
      L4400:
        e = exp (-t);
        x = 1. / t;
        y = t - 12.5;
        w1 = (((x * -.18784686463512 + .22991849164985) * x - .49893752514047)
              * x - 2.1916512131607e-5) * e + sqrt (x * .785398163397448);
        w2 = ((((((((((y * -6.22272689880615e-15 + 1.04126809657554e-13) * y
                      - 6.842418230913e-13) * y + 1.576841731919e-11) * y -
                    4.203948834175e-10) * y + 6.287255934781e-9) * y -
                  8.307159819228e-8) * y + 1.356478091922e-6) * y -
                2.08065576105639e-5) * y + 2.5239673033234e-4) * y -
              .00294484050194539) * y + .0601396183129168;
        w3 = ((((((((((((y * -4.1956914545948e-17 + 5.94344180261644e-16) * y
                        - 1.148797566469e-14) * y + 1.881303962576e-13) * y -
                      2.413554618391e-12) * y + 3.372127423047e-11) * y -
                    4.933988617784e-10) * y + 6.116545396281e-9) * y -
                  6.69965691739299e-8) * y + 7.52380085447161e-7) * y -
                8.08708393262321e-6) * y + 6.88603417296672e-5) * y -
              4.67067112993427e-4) * y + .00542313365864597;
        w4 = (((((((((((((y * 2.90401781000996e-18 - 4.63389683098251e-17) *
                         y + 6.274018198326e-16) * y -
                        8.936002188168e-15) * y + 1.194719074934e-13) * y -
                      1.45501321259466e-12) * y + 1.64090830181013e-11) * y -
                    1.71987745310181e-10) * y + 1.63738403295718e-9) * y -
                  1.39237504892842e-8) * y + 1.06527318142151e-7) * y -
                7.27634957230524e-7) * y + 4.12159381310339e-6) * y -
              1.74648169719173e-5) * y + 8.50290130067818e-5;
        wts[m] *= w1 - w2 - w3 - w4;
        wts[m + 1] *= w2;
        wts[m + 2] *= w3;
        wts[m + 3] *= w4;
        r1 = (((((((((((y * 4.94869622744119e-17 + 8.0356880573916e-16) * y -
                       5.599125915431e-15) * y - 1.378685560217e-13) * y +
                     7.006511663249e-13) * y + 1.30391406991118e-11) * y +
                   8.06987313467541e-11) * y - 5.20644072732933e-9) * y +
                 7.72794187755457e-8) * y - 1.61512612564194e-6) * y +
               4.15083811185831e-5) * y - 7.87855975560199e-4) * y +
            .0114189319050009;
        r2 = (((((((((((y * 4.89224285522336e-16 + 1.06390248099712e-14) * y
                       - 5.446260182933e-14) * y - 1.613630106295e-12) * y +
                     3.910179118937e-12) * y + 1.90712434258806e-10) * y +
                   8.78470199094761e-10) * y - 5.97332993206797e-8) * y +
                 9.25750831481589e-7) * y - 2.02362185197088e-5) * y +
               4.92341968336776e-4) * y - .00868438439874703) * y +
            .115825965127958;
        r3 = ((((((((((y * 6.12419396208408e-14 + 1.12328861406073e-13) * y -
                      9.051094103059e-12) * y - 4.781797525341e-11) * y +
                    1.660828868694e-9) * y + 4.499058798868e-10) * y -
                  2.519549641933e-7) * y + 4.97744404018e-6) * y -
                1.25858350034589e-4) * y + .00270279176970044) * y -
              .0399327850801083) * y + .433467200855434;
        r4 = (((((((((((y * 4.63414725924048e-14 - 4.72757262693062e-14) * y
                       - 1.001926833832e-11) * y + 6.074107718414e-11) * y +
                     1.576976911942e-9) * y - 2.01186401974027e-8) * y -
                   1.84530195217118e-7) * y + 5.02333087806827e-6) * y +
                 9.66961790843006e-6) * y - .00158522208889528) * y +
               .0280539673938339) * y - .278953904330072) * y +
            1.82835655238235;
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: 15 =< T < 20 */
      L4500:
        e = exp (-t);
        x = 1. / t;
        y = t - 17.5;
        w1 = ((x * .1962326414943 - .4969524146449) * x - 6.0156581186481e-5)
            * e + sqrt (x * .785398163397448);
        w2 = (((((((((((y * -1.865060577297e-16 + 1.16661114435809e-15) * y +
                       2.563712856363e-14) * y - 4.498350984631e-13) * y +
                     1.765194089338e-12) * y + 9.04483676345625e-12) * y +
                   4.98930345609785e-10) * y - 2.11964170928181e-8) * y +
                 3.98295476005614e-7) * y - 5.49390160829409e-6) * y +
               7.74065155353262e-5) * y - .00148201933009105) * y +
            .0497836392625268;
        w3 = (((((((((((y * -5.54451040921657e-17 + 2.68748367250999e-16) * y
                       + 1.349020069254e-14) * y - 2.507452792892e-13) * y +
                     1.944339743818e-12) * y - 1.29816917658823e-11) * y +
                   3.49977768819641e-10) * y - 8.67270669346398e-9) * y +
                 1.31381116840118e-7) * y - 1.36790720600822e-6) * y +
               1.1921069767316e-5) * y - 1.42181943986587e-4) * y +
            .00412615396191829;
        w4 = ((((((((((((y * -7.56882223582704e-19 + 7.53541779268175e-18) *
                        y - 1.157318032236e-16) * y +
                       2.411195002314e-15) * y - 3.601794386996e-14) * y +
                     4.082150659615e-13) * y - 4.289542980767e-12) * y +
                   5.086829642731e-11) * y - 6.35435561050807e-10) * y +
                 6.82309323251123e-9) * y - 5.63374555753167e-8) * y +
               3.57005361100431e-7) * y - 2.40050045173721e-6) * y +
            4.94171300536397e-5;
        wts[m] *= w1 - w2 - w3 - w4;
        wts[m + 1] *= w2;
        wts[m + 2] *= w3;
        wts[m + 3] *= w4;
        r1 = (((((((((((y * 4.36701759531398e-17 - 1.12860600219889e-16) * y
                       - 6.149849164164e-15) * y + 5.820231579541e-14) * y +
                     4.396602872143e-13) * y - 1.24330365320172e-11) * y +
                   6.71083474044549e-11) * y + 2.43865205376067e-10) * y +
                 1.67559587099969e-8) * y - 9.32738632357572e-7) * y +
               2.39030487004977e-5) * y - 4.68648206591515e-4) * y +
            .00834977776583956;
        r2 = (((((((((((y * 4.98913142288158e-16 - 2.60732537093612e-16) * y
                       - 7.775156445127e-14) * y + 5.766105220086e-13) * y +
                     6.4326967296e-12) * y - 1.39571683725792e-10) * y +
                   5.95451479522191e-10) * y + 2.42471442836205e-9) * y +
                 2.4748571014312e-7) * y - 1.14710398652091e-5) * y +
               2.71252453754519e-4) * y - .00496812745851408) * y +
            .082602060202678;
        r3 = (((((((((((y * 1.91498302509009e-15 + 1.48840394311115e-14) * y
                       - 4.316925145767e-13) * y + 1.186495793471e-12) * y +
                     4.615806713055e-11) * y - 5.54336148667141e-10) * y +
                   3.48789978951367e-10) * y - 2.79188977451042e-9) * y +
                 2.09563208958551e-6) * y - 6.76512715080324e-5) * y +
               .00132129867629062) * y - .0205062147771513) * y +
            .288068671894324;
        r4 = (((((((((((y * -5.43697691672942e-15 - 1.12483395714468e-13) * y
                       + 2.826607936174e-12) * y - 1.26673449328e-11) * y -
                     4.258722866437e-10) * y + 9.45486578503261e-9) * y -
                   5.86635622821309e-8) * y - 1.28835028104639e-6) * y +
                 4.41413815691885e-5) * y - 7.61738385590776e-4) * y +
               .0096609090298555) * y - .101410568057649) * y +
            .954714798156712;
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: 20 =< T < 35 */
      L4600:
        e = exp (-t);
        x = 1. / t;
        w1 = ((x * .1962326414943 - .4969524146449) * x - 6.0156581186481e-5)
            * e + sqrt (x * .785398163397448);
        w2 = ((((((t * 7.29841848989391e-4 - .0353899555749875) * t +
                  2.07797425718513) * t - 100.464709786287) * t +
                3152.06108877819) * t - 62705.4715090012) * t + (x *
                                                                 15472124.6264919
                                                                 -
                                                                 5260743.91316381)
              * x + 767135.400969617) * e + w1 * .234479815323517;
        w3 = ((((((t * 2.36392855180768e-4 - .00916785337967013) * t +
                  .462186525041313) * t - 19.694378600654) * t +
                499.169195295559) * t - 6214.1984584509) * t +
              ((x * 52144505.3212414 - 13411346.4389309) * x +
               1136732.98305631) * x - 2815.01182042707) * e +
            w1 * .0192704402415764;
        if (t <= 25.)
        {
            w4 = (((((((t * 2.33766206773151e-7 - 3.81542906607063e-5) * t +
                       .00351416601267) * t - .166538571864728) * t +
                     4.80006136831847) * t - 87.3165934223603) * t +
                   977.683627474638) * t + x * 16600.094511764 -
                  6144.79071209961) * e + w1 * 2.25229076750736e-4;
        }
        else
        {
            w4 = ((((((t * 5.74245945342286e-6 - 7.58735928102351e-5) * t +
                      2.35072857922892e-4) * t - .00378812134013125) * t +
                    .309871652785805) * t - 7.11108633061306) * t +
                  55.5297573149528) * e + w1 * 2.25229076750736e-4;
        }
        wts[m] *= w1 - w2 - w3 - w4;
        wts[m + 1] *= w2;
        wts[m + 2] *= w3;
        wts[m + 3] *= w4;
        r1 = ((((((t * -4.45711399441838e-5 + .00127267770241379) * t -
                  .236954961381262) * t + 15.4330657903756) * t -
                522.799159267808) * t + 10595.1216669313) * t + (x *
                                                                 -2511772.35556236
                                                                 +
                                                                 872975.373557709)
              * x - 129194.382386499) * e + .145303521503316 / (t -
                                                                .145303521503316);
        r2 = (((((t * -.0785617372254488 + 6.35653573484868) * t -
                 338.29693876399) * t + 12512.0495802096) * t -
               316847.570511637) * t + ((x * -1024274661.27427 +
                                         370104713.293016) * x -
                                        58711900.5093822) * x +
              5386142.11391604) * e + 1.33909728812636 / (t -
                                                          1.33909728812636);
        r3 = (((((t * -.237900485051067 + 18.4122184400896) * t -
                 1002.00731304146) * t + 37515.1841595736) * t -
               950626.66339013) * t + ((x * -2881390146.51985 +
                                        1066259150.44526) * x -
                                       172465289.687396) * x +
              16041939.0230055) * e + 3.92696350135829 / (t -
                                                          3.92696350135829);
        r4 = ((((((t * -6.00691586407385e-4 - .364479545338439) * t +
                  15.7496131755179) * t - 654.944248734901) * t +
                17083.0039597097) * t - 290517.939780207) * t +
              (x * 34905969.8304732 - 16494452.2586065) * x +
              2968179.40164703) * e + 8.58863568901199 / (t -
                                                          8.58863568901199);
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: 35 =< T < 53 */
      L4700:
        x = t * t;
        e = exp (-t) * x * x;
        w1 = sqrt (.785398163397448 / t);
        w2 = ((t * 6.16374517326469e-4 - .0126711744680092) * t +
              .0814504890732155) * e + w1 * .234479815323517;
        w3 = ((t * 2.0829496985723e-4 - .00377489954837361) * t +
              .0209857151617436) * e + w1 * .0192704402415764;
        w4 = ((t * 5.7663198200099e-6 - 7.8918728380489e-5) * t +
              3.28297971853126e-4) * e + w1 * 2.25229076750736e-4;
        wts[m] *= w1 - w2 - w3 - w4;
        wts[m + 1] *= w2;
        wts[m + 2] *= w3;
        wts[m + 3] *= w4;
        r1 = ((t * -4.075575259146e-5 - 6.88846864931685e-4) * t +
              .0174725309199384) * e + .145303521503316 / (t -
                                                           .145303521503316);
        r2 = ((t * -3.62569791162153e-4 - .00909231717268466) * t +
              .184336760556262) * e + 1.33909728812636 / (t -
                                                          1.33909728812636);
        r3 = ((t * -9.65842534508637e-4 - .0449822013469279) * t +
              .608784033347757) * e + 3.92696350135829 / (t -
                                                          3.92696350135829);
        r4 = ((t * -.00219135070169653 - .119108256987623) * t -
              .750238795695573) * e + 8.58863568901199 / (t -
                                                          8.58863568901199);
        rts[m] = r1 / (r1 + 1.);
        rts[m + 1] = r2 / (r2 + 1.);
        rts[m + 2] = r3 / (r3 + 1.);
        rts[m + 3] = r4 / (r4 + 1.);
        m += 4;
        goto L400;


/*             ...T-range: T >= 53 */
      L4800:
        w1 = sqrt (.785398163397448 / t);
        w2 = w1 * .234479815323517;
        w3 = w1 * .0192704402415764;
        w4 = w1 * 2.25229076750736e-4;
/*         R1 = R14 / (T - R14) */
/*         R2 = R24 / (T - R24) */
/*         R3 = R34 / (T - R34) */
/*         R3 = R44 / (T - R44) */
/*         RTS (M)   = R1 / (ONE + R1) */
/*         RTS (M+1) = R2 / (ONE + R2) */
/*         RTS (M+2) = R3 / (ONE + R3) */
/*         RTS (M+3) = R4 / (ONE + R4) */
        wts[m] *= w1 - w2 - w3 - w4;
        wts[m + 1] *= w2;
        wts[m + 2] *= w3;
        wts[m + 3] *= w4;
        rts[m] = .145303521503316 / t;
        rts[m + 1] = 1.33909728812636 / t;
        rts[m + 2] = 3.92696350135829 / t;
        rts[m + 3] = 8.58863568901199 / t;
        m += 4;
      L400:
        ;
    }
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
