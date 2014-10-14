
#pragma ident "@(#)__libmcr_tan.c 1.4 04/02/25 SMI"

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
 * correctly rounded tan
 * 1. tan(0) = 0 exact
 * 2. tan(inf/nan) = x-x
 * 3. tan(x) = x/(1-ulp/2) for |x|< 2**-1022
 * 4. tan(x) = 2**-300*(x*2**300 + 2**-1022) for |x|<= 2**-967
 * 5. tan(x) = x + 2**-1022 for  x < cbrt(6)*2**-27
 *    	Proof. Since x + e rounded to x for e < 1/2*ulp(x)
 *	we need only to find the maximum x such that x-tan(x) < 0.5*ulp_(x).
 *	Since tan(x)-x ~ x*x*x/3 for tiny x, it is easy to verify that if
 *		2**-27 <= x < cbrt(6)*2**-27, then
 *		ulp(x) = 2**(-27-52) = 2**-79, and
 *		tan(x)-x <~ x*x*x/3 < 2**-80 = 1/2*ulp(x).
 *	Since tan(x)-x is an increasing function, the inequality remain until
 *	x decrease to 2**-27-ulp, where ulp(x) = 2**-28. But by that time,
 *	tan(x)-x < x*x*x/3 = 2**(-81)/3 < 1/2 * 2**-80 = 1/2*ulp_(x).
 *	Thus, for x <= 1.353860343122586362e-08  (0x3E4D12ED 0x0AF1A27F)
 *	simply return x + 2**-1022.
 * 6. compute tan(x) by calling the almost correctly rounded __libmcr_mx_tan
 * 7. examine the result to see if it is nearly half-way case, if it is too
 *    close to the half-way case, call multiprecision routing __libmcr_mi_tan
 *    with increasing accuracy until it guarantee the correctly roundedness.
 */

#include "mcr.h"

#pragma weak tan = __libmcr_tan

static const double
twom1022 = 2.22507385850720138309e-308,	/* 0x00100000 00000000 = 2**-1022 */
two300 = 2.03703597633448608627e+90,	/* 0x52B00000 00000000 */
twom300 = 4.90909346529772655310e-91,	/* 0x2D300000 00000000 */
ohu = 9.99999999999999888978e-01,	/* 0x3fefffff ffffffff */
rndc = -1.02735137937997933477e+00;	/* 0xbff07007 ffff0000 */

double
__libmcr_tan(double x)
{
	int hx, ix, rnd = 0, nw, lx, ncrd;
	double z, er;
	float fw, ft;

	ix = hx = HIGH_WORD(x);
	ix = hx & 0x7fffffff;
	lx = LOW_WORD(x);
	if ((ix | lx) == 0)
		return (x);	/* tan(x) = x when x=0 */
	/* filter out exception argument */
	if (ix >= 0x7ff00000)
		return (x - x);	/* inf and nan */
	z = rndc;
	ft = (float) z;
	fw = ft * ft;
	/* if |x|<cbrt(6.0)*2**-27 */
	if (ix < 0x3E4D12ED || (ix == 0x3E4D12ED && lx < 0x0AF1A27F)) {
		z = twom1022;
		if (hx < 0)
			z = -z;
		if (ix >= 0x03800000)
			return (x + z);	/* if |x|>=2**-967 */
		else {		/* if |x|<2**-1022 */
			if (ix < 0x00100000)
				return (x / ohu);
			/* otherwise: 2**-1022 <= |x| < 2**-967 */
			else
				return ((x * two300 + z) * twom300);
		}
	}
	rnd = (*(int *) &fw) & 3;

	z = __libmcr_mx_tan(x, &er);
	ncrd = __libmcr_mx_check(&z, rnd, er);
	nw = 8;
	while (ncrd == 1) {
		z = __libmcr_mi_tan(x, nw, &ncrd, rnd);
		nw += 2;
	}
	return (z);
}
