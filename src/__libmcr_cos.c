
#pragma ident "@(#)__libmcr_cos.c 1.4 04/02/25 SMI"

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
 * correctly rounded cos(x)
 *
 * 1. cos(0) = 1 exact
 * 2. cos(inf/nan) = x-x
 * 3. cos(x) = 1-tiny  for x <= 1.05367121277235070127e-08 (3E46A09E 667F3BCC)
 *    	Proof. Since 1.0 - e rounded to 1.0 for e <= 2**-54 and 1.0-2**53 for
 * 	e > 2**-54, the rounded result of 1-2**-300 is equal to the rounded
 *	value of cos(x), provided 1-cos(x) <= 2**-54. The maximum value of |x|
 *	such that 1-cos(x) <= 2**-54 is sqrt(2**-53) chopped. This is obtained
 *	by noting 1-cos(x) = x*x/2 - O(x^4). Thus
 *		x <= sqrt(2**-53) chopped
 *		   = 1.05367121277235070127e-08 (3E46A09E 667F3BCC).
 * 4. compute cos(x) by calling the almost correctly rounded __libmcr_mx_cos
 * 5. examine the result to see if it is nearly half-way case, if it is too
 *    close to the half-way case, call multiprecision routing __libmcr_mi_cos
 *    with increasing accuracy until it guarantee the correctly roundedness.
 */

#include "mcr.h"

#pragma weak cos = __libmcr_cos

static const double
one = 1.0,			/* 0x3ff00000 00000000 */
twom300 = 4.90909346529772655310e-91,	/* 0x2D300000 00000000 */
rndc = -1.02735137937997933477e+00;	/* 0xbff07007 ffff0000 */

double
__libmcr_cos(double x)
{
	int hx, ix, rnd = 0, nw, lx, ncrd;
	double z, er;
	float fw, ft;

	hx = HIGH_WORD(x);
	ix = hx & 0x7fffffff;
	lx = LOW_WORD(x);
	if ((ix | lx) == 0)
		return (one);	/* cos(0) = 1 exact */
	z = rndc;

	/* filter out exception argument */
	if (ix >= 0x7ff00000)
		return (x - x);	/* inf and nan */
	ft = (float) z;
	er = 1.0;
	fw = ft * ft;

	/* if |x|<=sqrt(2**-53) chopped, return 1-tiny */
	if (ix <= 0x3e46a09e) {
		if (ix < 0x3e46a09e || ((unsigned) lx) <= (unsigned) 0x667f3bcc)
			return (er - twom300);	/* use variable er because */
						/* compiler may optimize the */
						/* constant expression */
						/* one-twom300 */
	}
	rnd = (*(int *) &fw) & 3;

	/* first call almost correctly rounded routine */
	z = __libmcr_mx_cos(x, &er);

	/*
	 * check for correctly roundedness, repeat calling special
	 * __libmcr_mi_cos if not
	 */
	ncrd = __libmcr_mx_check(&z, rnd, er);
	nw = 8;

	while (ncrd == 1) {
		z = __libmcr_mi_cos(x, nw, &ncrd, rnd);
		nw += 2;
	}
	return (z);
}
