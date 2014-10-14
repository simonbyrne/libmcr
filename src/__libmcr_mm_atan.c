
#pragma ident "@(#)__libmcr_mm_atan.c 1.5 04/02/25 SMI"

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
 * void __libmcr_mm_atan(mx,mc,nw)
 *
 * Kernel atan(mx) in mi format by argument reduction and Taylor series.
 * Assume mx finite with moderate size (less than one in absolute value).
 *
 * Method:
 *	1. Replace negative x by -x using atan(-mx)=-atan(mx).
 *	2. Base on
 *		atan(y)+atan(y) = atan(2y/(1-y*y))
 *	   we have, if we let x = 2y/(1-y*y), then  y = x/(1+sqrt(1+x*x)).
 *	   Thus for some small number t (say 1/8), we reduce x by:
 *	       w = x; i=0;
 *	       while (w>t) {
 *		    w = w/(1+sqrt(1+w*w)); i+=1;
 *	       }
 *
 *	   then atan(x) = 2**i * atan(w)
 *	   E.g., if x = inf, then
 *		y1 = inf/(1+sqrt(1+inf*inf)) = 1
 *		y2 = 1/(1+sqrt(2)) = sqrt(2)-1 = 0.4142135624...
 *		y3 = y2/(1+sqrt(1+y2*y2)) = 0.198912367...
 *		y4 = y3/(1+sqrt(1+y3*y3)) = 0.098491403...
 *
 *	3. Compute atan(w) by talyor serier
 *	   atan(w) = w - w^3/3 + w^5/5 - w^7/7 + ....
 *
 * Accuracy:
 *	Since one extra word is used, the number of reliable significant bits
 *	is 24*(nw-3).
 */

/*
 * threshold value for using Taylor Series: 0x200000	-->   1/8 0x100000
 * -->   1/16 0x80000		-->   1/32 ...
 */
#define	THRESHOLD 0x100000

#include <stdlib.h>

#include "mcr.h"


void
__libmcr_mm_atan(mx, mc, nw)
	int mx[], mc[], nw;
{
	int n, k, i, nwp, *mt, *mz, *mp, *mw, *mone, *mci;

	/* exceptions: atan(+-0)=+-0, atan(NaN)=NaN */

	/* if mx = 0 or > 0x7ff00000, return mx */
/* Redundant exception handling code: see __libmcr_atan */
#if 0
	if (mx[2] == 0 || mx[2] > 0x7ff00000) {
		if (mc != mx)
			for (i = 0; i < nw; i++)
				mc[i] = mx[i];
		return;
	}
#endif

	nwp = nw + 1;
	mt = (int *) calloc(nwp, sizeof (int));
	mz = (int *) calloc(nwp, sizeof (int));
	mp = (int *) calloc(nwp, sizeof (int));
	mw = (int *) calloc(nwp, sizeof (int));
	mci = (int *) calloc(nwp, sizeof (int));
	mone = (int *) calloc(nwp, sizeof (int));

	/* initialize mone = 1.0 */
	mone[0] = mw[0] = 1;
	mone[1] = 0;
	mone[2] = 1;
	for (i = 3; i < nwp; i++) {
		mone[i] = 0;
	}

	/* initialize n=0 and mw = mx (if mx=inf, then set mw=1 and n=1) */
	if (mx[2] == 0x7ff00000) {
		n = 1;
		mw[1] = 0;
		mw[2] = 1;
		for (i = 3; i < nwp; i++) {
			mw[i] = 0;
		}
	} else {
		n = 0;
		for (i = 1; i < nw; i++)
			mw[i] = mx[i];
		for (i = nw; i < nwp; i++)
			mw[i] = 0;
	}


	while (mw[1] >= 0 || (mw[1] == -1 && mw[2] >= THRESHOLD))
	{	/* w >= THRESHOLD */
		/* w = w/(1+sqrt(1+w*w)) */
		n += 1;
		__libmcr_mm_mul(mw, mw, mt, nwp);
		__libmcr_mm_add(mt, mone, mt, nwp);
		__libmcr_mm_sqrt(mt, mt, nwp);
		__libmcr_mm_add(mone, mt, mt, nwp);
		__libmcr_mm_div(mw, mt, mw, nwp);

	}

	/* compute  mc =  w^3/3 - w^5/5 + w^7/7 + .... */
	__libmcr_mm_mul(mw, mw, mz, nwp);	/* mz = mw*mw */
	__libmcr_mm_mul(mz, mw, mt, nwp);	/* mt = mw*mw*mw */
	__libmcr_mm_divi(mt, 3, mci, nwp);	/* mci= mw*mw*mw/3 */
	k = 3;
	mp[1] = mci[1];
	while (mw[1] - mp[1] < nwp - 1) {
		k = k + 2;
		__libmcr_mm_mul(mt, mz, mt, nwp);	/* mt = mw**k */
		if ((k & 2) != 0)
			i = k;
		else
			i = -k;
		/* mp = * (-1)**(k/2-1)*x**k/k */
		__libmcr_mm_divi(mt, i, mp, nwp);
		__libmcr_mm_add(mci, mp, mci, nwp);
	}
	/* atan(mw) = mw - mw^3/3 + mw^5/5 ... */
	__libmcr_mm_sub(mw, mci, mci, nw);

	/* atan(mx) =  2**n * atan(mw) */
	__libmcr_mm_scalbn(mci, nwp, n);
	mc[0] = mx[0];
	for (i = 1; i < nw; i++)
		mc[i] = mci[i];
	free(mt);
	free(mz);
	free(mp);
	free(mw);
	free(mci);
	free(mone);
}
