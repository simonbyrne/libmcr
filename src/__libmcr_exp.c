
#pragma ident "@(#)__libmcr_exp.c 1.5 04/02/25 SMI"

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

#ifdef DEBUG
#include <stdio.h>
#endif

/*
 * correctly rounded exp
 *
 * 1. exp(0) = 1 exact
 * 2. exp(nan) = nan
 * 3. exp(inf) = inf
 * 4. exp(x>o_threshold) = huge*huge
 * 5. exp(x<u_threshold) = tiny*tiny
 * 6. exp(x) = 1.0 + x  for  x < 2**-53
 * 7. if exp(x) is a subnormal number, use special addition to simulate
 *    the subnormal rounding.
 * 8. compute exp(x) by calling the almost correctly rounded __libmcr_mx_exp
 * 9. examine the result to see if it is nearly half-way case, if it is too
 *    close to the half-way case, call multiprecision routing
 *    __libmcr_mi_exp with increasing accuracy until it guarantee the
 *    correctly roundedness.
 */

#include "mcr.h"

// #pragma weak exp = __libmcr_exp

static const double
rndc = -1.02735137937997933477e+00,		/* 0xbff07007 ffff0000 */
o_threshold = 7.09782712893383973096e+02,	/* 0x40862E42, 0xFEFA39EF */
u_threshold = -7.45133219101941108420e+02,	/* 0x40874910, 0xD52D3051 */
huge = 1.0e300,
one = 1.0,
tiny = 2.22507385850720187716e-308;		/* 0x00100000, 0x00000001 */

double
__libmcr_exp(double x)
{
	int hx, ix, rnd = 0, nw, lx, ncrd, m, i, iz;
	double z, er, w, t;
	float fw, ft;

	ix = hx = HIGH_WORD(x);
	ix = hx & 0x7fffffff;
	lx = LOW_WORD(x);
	if ((ix | lx) == 0)
		return (one);	/* exp(0) = 1 */
	/* filter out non-finite argument */
	if (ix >= 0x40862E42) {	/* if |x|>=709.78... */
		if (ix >= 0x7ff00000) {
			if (((ix & 0xfffff) | lx) != 0)
				return (x + x);	/* NaN */
			else
			{
				/* exp(+-inf)={inf,0} */
				return ((hx >= 0) ? x : 0.0);
			}
		}
		if (x > o_threshold)
			return (huge * huge);	/* overflow */
		if (x < u_threshold)
			return (tiny * tiny);	/* underflow */
	}
	z = rndc;
	ft = (float) z;
	fw = ft * ft;

	/*
	 * if |x| < 2**-53, then the correctly rounded result of exp(x) is
	 * the same as the correctly rounded result of 1+x
	 */
	if (ix < 0x3ca00000)
		return (one + x);

	rnd = (*(int *) &fw) & 3;

	z = __libmcr_k_mx_exp(x, &er, &m);
	iz = HIGH_WORD(z);
	i = (iz >> 20) + m;
	if (i > 0) {
		ncrd = __libmcr_mx_check(&z, rnd, er);
		if (ncrd == 0)
			HIGH_WORD(z) = iz + (m << 20);
	} else {	/* subnormal output, becareful about double rounding */
		HIGH_WORD(w) = ((iz >> 20) - i + 1) << 20;
		LOW_WORD(w) = (int) (tiny * tiny); /* underflow inexact sig */
		t = w + z;	/* similate subnormal rounding by adding */
		er += z - (t - w);
		z = t;
		ncrd = __libmcr_mx_check(&z, rnd, er);
		/* mask off the additional leading bit from previous adding */
		if (ncrd == 0)
			HIGH_WORD(z) &= 0x000fffff;
	}
	nw = 8;
	while (ncrd == 1) {
		z = __libmcr_mi_exp(x, nw, &ncrd, rnd);
		nw += 2;
	}
	return (z);
}
