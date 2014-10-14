
#pragma ident "@(#)__libmcr_k_mm_exp.c 1.6 04/02/25 SMI"

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
 * double __libmcr_k_mm_exp(mx,mc,nw)
 *
 * Compute exp(mx) in mi format, where mx must be less than (1/2)*ln2
 * and in mi format.
 *
 * Method:
 *	1. Scale mx to mr by 2^-N, where N is choosen to make
 *		|mr|=|mx|/2^N < 1/2*ln2*(1/2)^-8 < 2^-1024
 *
 *	2. Compute exp(mr) by talyor serier
 *	   mc = 1+mr+mr^2/2!+mr^3/3!+mr^4/4!+...
 *
 *		Let mc = mx, mg = mr
 *		For i=2,...
 *			mg =  mg*mr/i
 *			if mc[2]-mg[2]>=nw stop  else
 *			mc = mc+mg
 *		Finally, mc = 1 + mc
 *
 *	3. Compute mc^(2^N) by squaring N times
 */

#include <stdlib.h>
#include "mcr.h"

#ifdef DEBUG
#include <stdio.h>
#define	H5(x) x[0], x[1], x[2], x[3], x[4]
#endif


void
__libmcr_k_mm_exp(mx, mc, nw)
	int mx[], mc[], nw;
{
	int n, k, j, i, *mt, *mr, *mz, *mci, nwp;
#ifdef DEBUG
	int mtail[4];
#endif

/* Redundant exception handling code: see __libmcr_exp */
#if 0
	if (mx[2] == 0) {	/* exp(0) = 1 */
		mc[0] = 1;
		mc[1] = 0;
		mc[2] = 1;
		for (i = 3; i < nw; i++)
			mc[i] = 0;
		return;
	}
	if (mx[2] >= 0x7ff00000) {
		if (mc[2] > 0x7ff00000) {	/* exp(NaN) = NaN */
			mc[0] = mx[0];
			mc[1] = 0;
			mc[2] = mx[2] | 0x80000;
			for (i = 3; i < nw; i++)
				mc[i] = 0;
		} else {
			mc[0] = 1;
			mc[1] = 0;
			for (i = 3; i < nw; i++)
				mc[i] = 0;
			if (mx[0] == 1) {	/* exp(inf)=inf */
				mc[2] = mx[2];
			} else { /* exp(-inf) = 0 */
				mc[2] = 0;
			}
		}
	}
#endif

	nwp = nw + 1;
	mt = (int *) calloc(nwp + 1, sizeof (int));
	mr = (int *) calloc(nwp + 1, sizeof (int));
	mz = (int *) calloc(nwp + 1, sizeof (int));
	mci = (int *) calloc(nwp + 1, sizeof (int));

	for (i = 0; i < nw; i++)
		mr[i] = mx[i];
	for (i = nw; i < nwp; i++)
		mr[i] = 0;

#ifdef DEBUG
	printf("mx = % 1.20e  %X %X %X %X %X\n",
		__libmcr_mi_mitod(mx, nw, mtail), H5(mx));
#endif
	/* scale mr by 2^k so than mr < 2^-1024 */
	n = mx[1];
	n = (n << 4) + (n << 3);	/* n = 24 * mx[1] */
	i = mx[2];
	while (i != 0) {
		i >>= 1;
		n += 1;
	}
	if (n < -10)
		k = 0;
	else {
		k = n + 11;
		__libmcr_mm_scalbn(mr, nwp, -k);
	}
#ifdef DEBUG
	if (k != 0) {
		printf("mr = 2^%d * mx = % 1.20e  %X %X %X %X %X\n",
			k, __libmcr_mi_mitod(mr, nwp, mtail), H5(mr));
	}
#endif

	/* Taylor series */
	for (i = 0; i < nwp; i++)
		mt[i] = mci[i] = mr[i];
	j = 2;
	while (mci[1] - mt[1] < nwp && mt[2] != 0) { /* mt = mt*r/j = r^j/j! */
		__libmcr_mm_mul(mt, mr, mt, nwp);
		__libmcr_mm_divi(mt, j, mt, nwp);
		__libmcr_mm_add(mci, mt, mci, nwp);
		j++;
	}
	for (i = 3; i < nw; i++)
		mz[i] = 0;
	mz[0] = 1;
	mz[1] = 0;
	mz[2] = 1;		/* mz = 1.0 */
	__libmcr_mm_add(mz, mci, mci, nwp);

#ifdef DEBUG
	printf("exp(mr) by Taylor Series:\n");
	printf("mci = % 1.20e  %X %X %X %X %X\n",
		__libmcr_mi_mitod(mci, nwp, mtail), H5(mci));
#endif

	/* square k times */
	for (i = 0; i < k; i++)
		__libmcr_mm_mul(mci, mci, mci, nwp);
	for (i = 0; i < nw; i++)
		mc[i] = mci[i];
#ifdef DEBUG
	printf("mci^(2^%d) = % 1.20e  %X %X %X %X %X\n",
		k, __libmcr_mi_mitod(mci, nwp, mtail), H5(mci));
#endif
	free(mt);
	free(mr);
	free(mz);
	free(mci);
}
