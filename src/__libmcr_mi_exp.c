
#pragma ident "@(#)__libmcr_mi_exp.c 1.6 04/02/25 SMI"

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
 * double __libmcr_mi_exp(double x,int nw,int *ncrd,rnd)
 *
 * Using multi-precision arithmetic, __libmcr_mi_exp return an almost
 * correctly rounded exp(x), with *ncrd=1 indicate if it it not guarantee
 * correctly rounded. Here nw is the number of integer words in the
 * multi-precision array. (Three words used for sign, exponent, and
 * leading word. Each addition integer represented an additional
 * 24-bits precision.)
 *
 * Method:
 *	1. Argument Reduction: given the input x, find r and integer k
 *	   and j such that
 *             	x = k*(ln2) + r,  |r| <= (1/2)*ln2 .
 *
 *	   Note that in order to compute r = x - k*ln2, we first examine
 *	   how close k*ln2 to a 53 double, and we find that when
 *	   n = 116877, n*ln2 is the closest one to a 53 bit floating
 *	   point number (only the significant parts are displayed)
 *           116877*ln2  1.B3C74F688A1382000033BCA50B1C3
 *           rounded to  1.B3C74F688A1382
 *                                       ^^^^^ 18 binary zeros
 *	   Thus, 3 extra words (3*24=72 bits) is needed to recover the
 *	   cancellation.
 *
 *	2. call kernel function __libmcr_k_mm_exp to compute exp(r) in mi format
 *
 *	3. rounded the mi format exponential to double and scale 2**k
 *
 * Special cases:
 *	exp(INF) is INF, exp(NaN) is NaN;
 *	exp(-INF)=  0;
 *	for finite argument, only exp(0)=1 is exact.
 *
 * Accuracy:
 *	according to an error analysis, the error is always somewhat less than
 *	2**-22 ulp (unit in the last place).
 *
 * Misc. info.
 *	For IEEE double
 *		if x >  7.09782712893383973096e+02 then exp(x) overflow
 *		if x < -7.45133219101941108420e+02 then exp(x) underflow
 *
 * Constants:
 * The hexadecimal values are the intended ones for the following constants.
 * The decimal values may be used, provided that the compiler will convert
 * from decimal to binary accurately enough to produce the hexadecimal values
 * shown.
 */

#include <stdlib.h>
#include "mcr.h"

static const double
half[] = {0.5, -0.5},
/* Redundant exception handling code: see __libmcr_exp */
#if 0
o_threshold = 7.09782712893383973096e+02,	/* 0x40862E42, 0xFEFA39EF */
u_threshold = -7.45133219101941108420e+02,	/* 0x40874910, 0xD52D3051 */
huge = 1.0e300,
one = 1.0,
twom1000 = 9.33263618503218878990e-302,
#endif
invln2 = 1.44269504088896338700e+00;	/* 0x3ff71547, 0x652b82fe */

double
__libmcr_mi_exp(x, nw, ncrd, rnd)
	double x;
	int nw, *ncrd, rnd;
{
	/*
	 * rnd: rounding direction: 0 --  nearest, 1 - chopped, 2 - up, 3 -
	 * down
	 */

	int xsb, n, *mc, *mr, *mx, *my, *ln2, nwp;
	int hx;
	double z;

	hx = HIGH_WORD(x);	/* high word of x */
	xsb = (hx >> 31) & 1;	/* sign bit of x */
	hx = hx & 0x7fffffff;	/* high word of |x| */

	/* filter out non-finite argument */
	*ncrd = 0;
/* Redundant exception handling code: see __libmcr_exp */
#if 0
	if (hx >= 0x40862E42) {	/* if |x|>=709.78... */
		if (hx >= 0x7ff00000) {
			if (((hx & 0xfffff) | LOW_WORD(x)) != 0)
				return ((x + x)); /* NaN */
			else
			{
				/* exp(+-inf)={inf,0} */
				return (((xsb == 0) ? x : 0.0));
			}
		}
		if (x > o_threshold)
			return (huge * huge);	/* overflow */
		if (x < u_threshold)
			return (twom1000 * twom1000);	/* underflow */
	}
	if (hx <= 0x3c900000) {	/* |x|<= 2**-54 */
		return (one + x);
	}
#endif
	n = x * invln2 + half[xsb];

	nwp = nw + 3;
	mx = (int *) calloc(nwp, sizeof (int));
	my = (int *) calloc(nwp, sizeof (int));
	mr = (int *) calloc(nwp, sizeof (int));
	ln2 = (int *) calloc(nwp, sizeof (int));
	__libmcr_mm_ln2(ln2, nwp);

	mc = (int *) calloc(nw, sizeof (int));

	/* compute  r = x - n*ln2 */
	__libmcr_mi_dtomi(x, mx, nwp);

	__libmcr_mm_muli(ln2, n, my, nwp);
	__libmcr_mm_sub(mx, my, mr, nwp);

	/* compute mc = exp(mr) */
	__libmcr_k_mm_exp(mr, mc, nw);

	/* scale mr by 2**n */
	__libmcr_mm_scalbn(mc, nw, n);

	/* final rounding */
	z = __libmcr_mi_final(mc, nw, ncrd, rnd);
	free(mc);
	free(mr);
	free(mx);
	free(my);
	free(ln2);
	return (z);
}
