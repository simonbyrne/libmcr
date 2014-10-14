
#pragma ident "@(#)__libmcr_mx_tan.c 1.5 04/02/25 SMI"

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
 * double __libmcr_mx_tan((double)x,(double*)e)
 *
 * Method:
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

#include "mcr.h"

double
__libmcr_mx_tan(x, e)
	double x, *e;
{
	double r, y[2], v, r_h, v_h, v_l;
	int ix, n;

	/* High word of x. */
	ix = HIGH_WORD(x);

	/* |x| ~< pi/4 */
	ix &= 0x7fffffff;
	if (ix <= 0x3fe921fb) {
		r = __libmcr_k_mx_tan(x, 0.0, e);
	}
/* Redundant exception handling code: see __libmcr_tan */
#if 0
	/* trig(Inf or NaN) is NaN */
	else if (ix >= 0x7ff00000)
		return (x - x);
#endif

	/* argument reduction needed */
	else {
		n = __libmcr_k_mx_rem_pio2(x, y);
		r = __libmcr_k_mx_tan(y[0], y[1], e);
		if ((n & 1) == 1) {	/* return -1/r  */
			v = -1.0 / r;
			r_h = (double) ((float) r);
			v_h = (double) ((float) v);
			v_l = 1.0 + r_h * v_h;
			v_l += v_h * (*e - (r_h - r));
			v_l *= v;
			r = v_h + v_l;
			*e = v_l - (r - v_h);
		}
	}
	return (r);
}
