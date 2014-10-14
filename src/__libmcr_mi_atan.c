
#pragma ident "@(#)__libmcr_mi_atan.c 1.6 04/02/25 SMI"

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
 * double __libmcr_mi_atan(double x,int nw,int *ncrd,rnd)
 *
 * Using multi-precision arithmetic, __libmcr_mi_atan return an almost
 * correctly rounded atan(x), with *ncrd=1 to indicate it it not guaranteed
 * correctly rounded. Here nw is the number of integer words in the
 * multi-precision array. (Three words used for sign, exponent, and
 * leading word. Each addition integer represented an additional
 * 24-bits precision.)
 *
 * Method: (the algorithm is based on atan(x)=atan(y)+atan((x-y)/(1+x*y))).
 *	or  atan(x) + atan(y) = atan((x+y)/(1-x*y))
 *
 *	Assume __libmcr_mm_atan(x) return atan(x) for x<sqrt(2)-1.
 *
 *	1. reduce x to |x| by atan(-x) = atan(x) if necessary.
 *	2. If x <= sqrt(2)-1 return __libmcr_mm_atan(x)
 *	3. If x >= sqrt(2)+1 return atan(inf) - __libmcr_mm_atan(1/x)
 *	4. Otherwise:
 *		atan(x) = atan(1.0) + __libmcr_mm_atan((x-1)/(x+1))
 *
 *	5. determine whether the error is close to 0.5 ulp by applying
 *	   __libmcr_mi_final()
 *
 */

#include <stdlib.h>
#include "mcr.h"

/* Redundant exception handling code: see __libmcr_atan */
#if 0
static const double
pio2h = 1.570796326794896558e+00,
pio2l = 6.123233995736765886e-17;
#endif

double
__libmcr_mi_atan(x, nw, ncrd, rnd)
	double x;
	int nw, *ncrd, rnd;
{
	/*
	 * rnd: rounding direction: 0 --  nearest, 1 - chopped, 2 - up, 3 -
	 * down
	 */
	double z;
	int *mc, *mr, *mx, *my, *pio2, *pio4;
	int hx, ix, i;
/* Redundant exception handling code: see __libmcr_atan */
#if 0
	int lx;
#endif

	hx = HIGH_WORD(x);	/* high word of x */
	ix = hx & 0x7fffffff;	/* high word of |x| */

	/* filter out special argument */
	*ncrd = 0;
/* Redundant exception handling code: see __libmcr_atan */
#if 0
	lx = LOW_WORD(x);	/* low word of x */
	if ((ix | lx) == 0)
		return (x);	/* 0 */
	if (ix >= 0x7ff00000) {
		if (((ix - 0x7ff00000) | lx) != 0)
			return (x + x);	/* NaN */
		/* inf: return pi with sign */
		z = -pio2l;
		if (hx < 0)
			return (z - pio2h);	/* atan(-inf) = -pi/2 */
		else
			return (pio2h - z);	/* atan(inf) = pi/2 */
	}
#endif

	mc = (int *) calloc(nw, sizeof (int));
	mr = (int *) calloc(nw, sizeof (int));
	mx = (int *) calloc(nw, sizeof (int));
	my = (int *) calloc(nw, sizeof (int));
	pio2 = (int *) calloc(nw, sizeof (int));
	pio4 = (int *) calloc(nw, sizeof (int));
	__libmcr_mm_pio2(pio2, nw);
	__libmcr_mm_pio4(pio4, nw);

	__libmcr_mi_dtomi(x, mx, nw);
	mx[0] = 1.0;
	if (ix <= 0x3FDA8279)	/* |x| < 0.4142... */
		__libmcr_mm_atan(mx, mc, nw);
	else {
		my[0] = 1;
		my[1] = 0;
		my[2] = 1.0;
		for (i = 3; i < nw; i++)
			my[i] = 0;	/* my= +1 */
		if (ix >= 0x4003504F) {	/* |x| > sqrt(2)+1= 2.4142 */
			__libmcr_mm_div(my, mx, mr, nw);
			__libmcr_mm_atan(mr, mr, nw);
			/* mc = atan(inf) - atan(1/x) */
			__libmcr_mm_sub(pio2, mr, mc, nw);
		} else {	/* sqrt(2)+1 > |x| > sqrt(2)-1 =.4142 */
			/*
			 * atan(x) = atan(1.0) +
			 * __libmcr_mm_atan((x-1)/(x+1))
			 */
			__libmcr_mm_sub(mx, my, mr, nw);	/* mr = x-1 */
			__libmcr_mm_add(mx, my, my, nw);	/* my = x+1 */
			__libmcr_mm_div(mr, my, mr, nw);
			__libmcr_mm_atan(mr, mr, nw);
			__libmcr_mm_add(pio4, mr, mc, nw);
		}
	}
	if (hx < 0)
		mc[0] = -mc[0];

	/* final rounding */
	z = __libmcr_mi_final(mc, nw, ncrd, rnd);
	free(mc);
	free(mr);
	free(mx);
	free(my);
	free(pio2);
	free(pio4);
	return (z);
}
