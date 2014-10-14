
#pragma ident "@(#)__libmcr_mi_log.c 1.5 04/02/25 SMI"

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
 * double __libmcr_mi_log(double x,int nw,int *ncrd,rnd)
 *
 * Using multi-precision arithmetic, __libmcr_mi_log return an almost
 * correctly rounded log(x), with *ncrd=1 indicate if it it not guarantee
 * correctly rounded. Here nw is the number of integer words in the
 * multi-precision array. (Three words used for sign, exponent, and
 * leading word. Each addition integer represented an additional
 * 24-bits precision.)
 *
 * Method:
 *	1. Argument Reduction: given the input x, find n such that
 *             	x = 2**n * r,  1/sqrt2 <= r <= sqrt2
 *
 *	2. call __libmcr_k_mi_log to compute log(r)
 *
 *	3. Log(x) = n*ln2 + log(r)
 *
 * Special cases:
 *	log(INF) is INF, log(NaN) is NaN;
 *	log(0) is -INF (with divide-by-zero signal);
 *	log(x<0)=  NaN (with invalid signal);
 *	for finite argument, only log(1)=0 is exact.
 *
 * Accuracy:
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
/* Redundant exception handling code: see __libmcr_log */
#if 0
one = 1.0,
zero = 0.0,
#endif
two54 = 18014398509481984.0;

double
__libmcr_mi_log(x, nw, ncrd, rnd)
	double x;
	int nw, *ncrd, rnd;
{
	/*
	 * rnd: rounding direction: 0 --  nearest, 1 - chopped, 2 - up, 3 -
	 * down
	 */

	int k, *mc, *mt, *ln2;
	int hx, na;
	double z;
/* Redundant exception handling code: see __libmcr_log */
#if 0
	int i;
#endif

	hx = HIGH_WORD(x);	/* high word of x */
	na = 0;

/* Redundant exception handling code: see __libmcr_log */
#if 0
	/* filter out non-finite argument */
	if ((hx + 0x00208000) < 0x00308000) { /* subnorm,0,neg,inf,nan,huge */
		i = ((unsigned) (hx << 1)) >> 1;
		if (hx >= 0x7ff00000)
			return (x + x);	/* x is inf/nan */
		if ((i | LOW_WORD(x)) == 0)
			return (-one / (x * x));	/* log(+-0) is -inf */
		if (hx < 0)
			return ((x - x) / zero);	/* log(-#) is NaN */
#endif
		if (hx < 0x00100000) {	/* subnormal */
			x *= two54;
			hx = HIGH_WORD(x);
			na = -54;
		}
/* Redundant exception handling code: see __libmcr_log */
#if 0
	}
#endif
	k = ((hx + 0x95f61) >> 20) - 0x3ff;
	na += k;
	hx -= (k << 20);	/* normalize x to [1/sqrt2, sqrt2] */
	HIGH_WORD(x) = hx;

	mc = (int *) calloc(nw, sizeof (int));
	mt = (int *) calloc(nw, sizeof (int));
	ln2 = (int *) calloc(nw, sizeof (int));
	__libmcr_mm_ln2(ln2, nw);

	/* call __libmcr_k_mi_log */
	__libmcr_k_mi_log(x, mc, nw);

	/* 3. Log(x) = n*ln2 + log(r) */
	__libmcr_mm_muli(ln2, na, mt, nw);
	__libmcr_mm_add(mt, mc, mc, nw);
	z = __libmcr_mi_final(mc, nw, ncrd, rnd);
	free(mt);
	free(mc);
	free(ln2);
	return (z);
}
