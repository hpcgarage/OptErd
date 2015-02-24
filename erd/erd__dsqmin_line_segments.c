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

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__DSQMIN_LINE_SEGMENTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This function returns the square of the minimum 3D */
/*                distance that exists between two line segments P and Q */
/*                in space. */
/*                The line segment P is located between points P0 and */
/*                P1 and likewise for line segment Q between points */
/*                Q0 and Q1. The routine also handles the situations */
/*                when the line segments P and/or Q have zero lenghts, */
/*                i.e. they are points in space. */
/*                        \ */
/*                         \ */
/*                          P0 */
/*                           \ <-- segment P, dir vector U = P1 - P0 */
/*                            \ */
/*                             P1 */
/*                              \ */
/*                               \          / */
/*                                \        / <-- line L2 */
/*                                 P(sc)  / */
/*                                 |\    / */
/*                                 | \  / */
/*             W = P0 - Q0         |  \/ */
/*      C(sc,tc) = P(sc) - Q(tc)   |  /\ */
/*                                 | /  \ */
/*                                 |/    \ */
/*                                 Q(tc)  \ */
/*                                /        \ */
/*                               /          \ <-- line L1 */
/*                              /            \ */
/*                             Q1             \ */
/*                            / */
/*                           / <-- segment Q, dir vector V = Q1 - Q0 */
/*                          Q0 */
/*                         / */
/*                The two lines L1 and L2 are the extension lines on */
/*                which the two segments P and Q sit and are in given */
/*                in parametric form using the direction vectors U */
/*                and V as: */
/*                                 L1 = P0 + sU */
/*                                 L2 = Q0 + tV */
/*                with s and t parameters. C(sc,tc) is the vector */
/*                representing the closest approach between these two */
/*                lines occuring at s = sc and t = tc and is uniquely */
/*                perpendicular to the line direction vectors U and V: */
/*                                 U dot C(sc,tc) = 0 */
/*                                 V dot C(sc,tc) = 0 */
/*                Using C(sc,tc) = P(sc) - Q(tc) = P0 + sc*U - Q0 - tc*V */
/*                = W + sc*U - tc*V for these two conditions we obtain */
/*                two simultaneous linear equations in sc and tc: */
/*                     (U dot U)*sc - (U dot V)*tc = - (U dot W) */
/*                     (V dot U)*sc - (V dot V)*tc = - (V dot W) */
/*                from which we obtain: */
/*                          sc = (b*e - c*d) / (a*c - b*b) */
/*                          tc = (a*e - b*d) / (a*c - b*b) */
/*                where: */
/*                                 a = U dot U */
/*                                 b = U dot V */
/*                                 c = V dot V */
/*                                 d = U dot W */
/*                                 e = V dot W */
/*                If the denominator (a*c - b*b) is equal to zero, then */
/*                the lines L1 and L2 are parallel and the distance */
/*                between them is constant. We can solve for this */
/*                parallel distance of separation by fixing the value of */
/*                one parameter and using either equation to solve for */
/*                the other. If we select sc = 0, then we get tc = d/b */
/*                = e/c. The situation when one or both of the line */
/*                segments shrink to points is identified as a and/or c */
/*                being equal to zero. */

/*                The closest distance between segments may not be the */
/*                same as the closest distance between their extended */
/*                lines. The closest points on the extended infinite line */
/*                may be outside the range of the segments. Using the */
/*                parametrized line form for L1 = P0 + sU, we see that */
/*                the segment P is represented by all points on L1 for */
/*                which 0 =< s =< 1. Likewise for the segment Q with */
/*                0 =< t =< 1. */
/*                Minimizing the distance between the two line segments */
/*                P and Q means to minimize the scalar product: */
/*                  C dot C = (W0 + s*U - t*V) dot (W0 + s*U - t*V) */
/*                as a function of s and t over the restricted values */
/*                0 =< s =< 1 and 0 =< t =< 1 defining the segments. */
/*                  Input: */
/*                     XP0,YP0,ZP0 =  x,y,z coordinates of first point */
/*                                    defining line segment P */
/*                     XP1,YP1,ZP1 =  x,y,z coordinates of second point */
/*                                    defining line segment P */
/*                     XQ0,YQ0,ZQ0 =  x,y,z coordinates of first point */
/*                                    defining line segment Q */
/*                     XQ1,YQ1,ZQ1 =  x,y,z coordinates of second point */
/*                                    defining line segment Q */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD double erd__dsqmin_line_segments (double xp0, double yp0,
                                  double zp0, double xp1,
                                  double yp1, double zp1,
                                  double xq0, double yq0,
                                  double zq0, double xq1,
                                  double yq1, double zq1)
{
    /* ...set direction vectors U and V and reference vector W. */
    const double xu = xp1 - xp0;
    const double yu = yp1 - yp0;
    const double zu = zp1 - zp0;
    const double xv = xq1 - xq0;
    const double yv = yq1 - yq0;
    const double zv = zq1 - zq0;
    const double xw = xp0 - xq0;
    const double yw = yp0 - yq0;
    const double zw = zp0 - zq0;
    const double udotu = xu * xu + yu * yu + zu * zu;
    const double vdotv = xv * xv + yv * yv + zv * zv;
    double xc, yc, zc;
    if (udotu < 1.0e-12 && vdotv < 1.0e-12) {
        /* ...both line segments are points. */
        return xw * xw + yw * yw + zw * zw;
    } else if (udotu < 1.0e-12) {
        /* ...only P line segment is a point. */
        const double vdotw = xv * xw + yv * yw + zv * zw;
        if (vdotw < 0.0) {
            return xw * xw + yw * yw + zw * zw;
        } else if (vdotw > vdotv) {
            xc = xw - xv;
            yc = yw - yv;
            zc = zw - zv;
        } else {
            const double tc = vdotw / vdotv;
            xc = xw - tc * xv;
            yc = yw - tc * yv;
            zc = zw - tc * zv;
        }
    } else if (vdotv < 1.0e-12) {
        /* ...only Q line segment is a point. */
        const double udotw = -(xu * xw + yu * yw + zu * zw);
        if (udotw < 0.0) {
            return xw * xw + yw * yw + zw * zw;
        } else if (udotw > udotu) {
            xc = xw + xu;
            yc = yw + yu;
            zc = zw + zu;
        } else {
            const double sc = udotw / udotu;
            xc = xw + sc * xu;
            yc = yw + sc * yu;
            zc = zw + sc * zu;
        }
    } else {
        /* ...both P and Q are line segments. */
        const double udotv = xu * xv + yu * yv + zu * zv;
        const double udotw = xu * xw + yu * yw + zu * zw;
        const double vdotw = xv * xw + yv * yw + zv * zw;
        const double denom = udotu * vdotv - udotv * udotv;
        double sn, sd, tn, td;
        if (denom < 1.0e-12) {
            sn = 0.0;
            sd = 1.0;
            tn = vdotw;
            td = vdotv;
        } else {
            sn = udotv * vdotw - vdotv * udotw;
            sd = denom;
            if (sn < 0.0) {
                sn = 0.0;
                tn = vdotw;
                td = denom;
            } else if (sn > denom) {
                sn = denom;
                sd = denom;
                tn = vdotw + udotv;
                td = vdotv;
            } else {
                tn = udotu * vdotw - udotv * udotw;
                td = denom;
            }
        }
        if (tn < 0.0) {
            tn = 0.0;
            sn = -udotw;
            if (sn < 0.0) {
                sn = 0.0;
            } else if (sn > udotu) {
                sn = sd;
            } else {
                sn = -udotw;
                sd = udotu;
            }
        } else if (tn > td) {
            tn = td;
            sn = udotv - udotw;
            if (sn < 0.0) {
                sn = 0.0;
            } else if (sn > udotu) {
                sn = sd;
            } else {
                sd = udotu;
            }
        }
        const double sc = sn / sd;
        const double tc = tn / td;
        xc = xw + sc * xu - tc * xv;
        yc = yw + sc * yu - tc * yv;
        zc = zw + sc * zu - tc * zv;
    }
    return xc * xc + yc * yc + zc * zc;
}
