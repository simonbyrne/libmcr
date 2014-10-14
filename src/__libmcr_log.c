
#pragma ident "@(#)__libmcr_log.c 1.4 04/02/25 SMI"

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
 * correctly rounded log
 *
 * 1. log(1) = 0, log(inf)=inf exact
 * 2. log(nan) = nan
 *    log(-ve) = NaN (invalid exception)
 *    log(0) = -inf (divide-by-zero exception)
 * 4. compute log(x) by calling the almost correctly rounded __libmcr_mx_log
 * 5. examine the result to see if it is nearly half-way case, if it is too
 *    close to the half-way case, call multiprecision routing __libmcr_mi_log
 *    with increasing accuracy until it guarantee the correctly roundedness.
 */

#include "mcr.h"

//#pragma weak log = __libmcr_log

static const double
one = 1.0,
rndc = -1.02735137937997933477e+00;	/* 0xbff07007 ffff0000 */

double
__libmcr_log(double x)
{
	int hx, ix, rnd, nw, lx, ncrd;
	double z, er;
	float fw, ft;

	ix = hx = HIGH_WORD(x);
	ix = hx & 0x7fffffff;
	lx = LOW_WORD(x);
	if (((hx - 0x3ff00000) | lx) == 0)
		return (0.0);	/* log(1) = 0 */
	/* filter out exception argument */
	if ((hx + 0x00208000) < 0x00308000) {
		/* subnormal,0,negative,inf,na,n,huge */
		if (ix >= 0x7ff00000)
			return (x + x);	/* x is inf/nan */
		if ((ix | lx) == 0)
			return (-one / (x * x));	/* log(+-0) is -inf */
		if (hx < 0)
			return ((x - x) / (x - x));	/* log(-#) is NaN */
	}
	z = rndc;
	ft = (float) z;
	fw = ft * ft;

	rnd = (*(int *) &fw) & 3;

	z = __libmcr_mx_log(x, &er);
	ncrd = __libmcr_mx_check(&z, rnd, er);
	nw = 8;
	while (ncrd == 1) {
		z = __libmcr_mi_log(x, nw, &ncrd, rnd);
#ifdef DEBUG
		printf("nw = %d\n", nw);
		printf("z = %08X %08X  ncrd = %d rnd=%d\n",
			HIGH_WORD(z), LOW_WORD(z), ncrd, rnd);
		fflush(stdout);
#endif
		nw += 2;
	}
	return (z);
}
