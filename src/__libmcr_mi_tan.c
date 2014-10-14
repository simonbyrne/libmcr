
#pragma ident "@(#)__libmcr_mi_tan.c 1.5 04/02/25 SMI"

/* ************************************************************ */
/*
 * Copyright (c) 2002 Sun Microsystems, Inc. All  Rights Reserved.
 *
 * Redistribution  and use in  source  and binary  forms,  with or
 * without modification, are permitted provided that the following
 * conditions are met:
 *
 *  -Redistribution of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 *  -Redistribution  in  binary  form  must  reproduce  the  above
 *   copyright notice,  this list of conditions and  the following
 *   disclaimer  in  the  documentation  and/or  other   materials
 *   provided with the distribution.
 *
 * Neither  the  name  of  Sun Microsystems, Inc.  or the names of
 * contributors may be used to endorse or promote products derived
 * from this  software without  specific prior written permission.
 *
 * This  software  is provided  "AS IS," without a warranty of any
 * kind.  ALL EXPRESS  OR IMPLIED  CONDITIONS, REPRESENTATIONS AND
 * WARRANTIES,  INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
 * FITNESS  FOR  A  PARTICULAR  PURPOSE  OR  NON-INFRINGEMENT, ARE
 * HEREBY  EXCLUDED.    SUN MICROSYSTEMS, INC.  ("SUN")   AND  ITS
 * LICENSORS  SHALL  NOT  BE  LIABLE  FOR  ANY DAMAGES SUFFERED BY
 * LICENSEE  AS A  RESULT OF USING, MODIFYING OR DISTRIBUTING THIS
 * SOFTWARE  OR  ITS  DERIVATIVES.   IN  NO  EVENT WILL SUN OR ITS
 * LICENSORS BE LIABLE  FOR  ANY LOST REVENUE,  PROFIT OR DATA, OR
 * FOR DIRECT,  INDIRECT,  SPECIAL,  CONSEQUENTIAL,  INCIDENTAL OR
 * PUNITIVE DAMAGES,  HOWEVER CAUSED AND  REGARDLESS OF THE THEORY
 * OF LIABILITY,  ARISING  OUT  OF  THE USE OF OR INABILITY TO USE
 * THIS SOFTWARE,  EVEN IF SUN HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES.
 *
 * You acknowledge that this software is not designed, licensed or
 * intended  for  use  in  the  design, construction, operation or
 * maintenance of any nuclear facility.
 *								*/
/* ************************************************************ */

/*
 * double __libmcr_mi_tan(double x,int nw,int *ncrd,int rnd)
 *
 * Using multi-precision arithmetic, __libmcr_mi_tan return an almost
 * correctly rounded tan(x), with *ncrd=1 indicate if it it not guarantee
 * correctly rounded. Here nw is the number of integer words in the
 * multi-precision array. (Three words used for sign, exponent, and
 * leading word. Each addition integer represented an additional
 * 24-bits precision.)
 *
 * kernel function:
 *      __libmcr_k_mm_sin         ... sine function on [-pi/4,pi/4]
 *      __libmcr_k_mm_cos         ... cose function on [-pi/4,pi/4]
 *      __libmcr_k_mi_rem_pio2    ... argument reduction routine
 *
 * Method.
 *      Let S and C denote the sin and cos respectively on [-PI/4, +PI/4].
 *      1. Assume the argument x is reduced to y1+y2 = x-k*pi/2 in
 *         [-pi/2 , +pi/2], and let n = k mod 4.
 *      2. Let S=S(y1+y2), C=C(y1+y2). Depending on n, we have
 *
 *          n        sin(x)      cos(x)        tan(x)
 *     ----------------------------------------------------------
 *          0          S           C             S/C
 *          1          C          -S            -C/S
 *          2         -S          -C             S/C
 *          3         -C           S            -C/S
 *     ----------------------------------------------------------
 *
 * Special cases:
 *      Let trig be any of sin, cos, or tan.
 *      trig(+-INF)  is NaN, with signals;
 *      trig(NaN)    is that NaN;
 */

#include <stdlib.h>
#include "mcr.h"

double
__libmcr_mi_tan(x, nw, ncrd, rnd)
	double x;
	int nw, *ncrd, rnd;
{
	/*
	 * rnd: rounding direction: 0 --  nearest, 1 - chopped, 2 - up, 3 -
	 * down
	 */

	int n, *ms, *mx, *mc, *mt;
	int ix;
	double z;

	ms = (int *) calloc(nw, sizeof (int));
	mx = (int *) calloc(nw, sizeof (int));
	mc = (int *) calloc(nw, sizeof (int));
	mt = (int *) calloc(nw, sizeof (int));

	/* High word of x. */
	ix = HIGH_WORD(x);

	/* |x| ~< pi/4 */
	ix &= 0x7fffffff;
	if (ix <= 0x3fe921fb) {
		__libmcr_mi_dtomi(x, mx, nw);
		__libmcr_k_mm_sin(mx, ms, nw);
		__libmcr_k_mm_cos(mx, mc, nw);
		__libmcr_mm_div(ms, mc, mt, nw);
	}
	/* trig(Inf or NaN) is NaN */
/* Redundant exception handling code: see __libmcr_tan */
#if 0
	else if (ix >= 0x7ff00000) {
		mt[0] = 1;
		mt[1] = 0;
		mt[2] = 0x7fffffff;
		for (n = 3; n < nw; n++)
			mt[n] = 0;
	}
#endif

	/* argument reduction needed */
	else {
		n = __libmcr_k_mi_rem_pio2(x, mx, nw);
		__libmcr_k_mm_sin(mx, ms, nw);
		__libmcr_k_mm_cos(mx, mc, nw);
		switch (n & 1) {
		case 0:
			__libmcr_mm_div(ms, mc, mt, nw);
			break;
		case 1:
			__libmcr_mm_div(mc, ms, mt, nw);
			mt[0] = -mt[0];
			break;
		}
	}

	/* final rounding */
	z = __libmcr_mi_final(mt, nw, ncrd, rnd);
	free(ms);
	free(mx);
	free(mc);
	free(mt);
	return (z);
}
