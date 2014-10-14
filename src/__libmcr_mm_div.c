
#pragma ident "@(#)__libmcr_mm_div.c 1.4 04/02/25 SMI"

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
 * mm_div compute mi number mx/my, return in mi number mq, nw the
 * number of words
 *
 * myi = my
 * mq  = mx
 * mri ~ 1/my (obtained by double precision arithmetic)
 * loop until mri=1+tiny, where tiny is less than 1/2 of the mantissa
 *	    myi = myi*mri
 *	    mq  = mq *mri
 *	    mri = 2-myi
 * Finally return q*r (the last update);
 *
 * Analysis:
 *	r0 ~ 1/y, q0 = x,  y0 = y
 * For i = 0, 1, 2,...
 *	y(i+1) = yi*ri
 *	q(i+1) = qi*ri
 *	r(i+1) = 2-y(i+1)
 *
 * If r0 = 1/y*(1+E0), then
 *	y1 = y0*r0 = (1+E0)
 * 	q1 = q0*r0 = x*r0 = x/y*(1+E0);
 *	r1 = 2 - y1 = 1 - E0.
 *
 *	y2 = y1*r1 = (1+E0)*(1-E0) = 1-E0^2
 *	q2 = q1*r1 = x/y*(1+E0)*(1-E0) = x/y*(1-E0^2)
 *	r2 = 2 - y2 = 1 + E0^2
 *
 *	...
 *
 * Accuracy:
 *	Since there are two extra words, the reliable is about 24*(nw-3).
 */
#include <stdlib.h>
#include "mcr.h"

#ifdef DEBUG
#include <stdio.h>
#endif

static const double
two24 = 16777216,
twom24 = 1.0 / 16777216.0;

void
__libmcr_mm_div(mx, my, mz, nw)
	int mx[], my[], mz[], nw;
{
	int i, k, j, mx2, my2, *myi, *mri, *mi_two, *mq, nwp;
#ifdef DEBUG
	int mtail[4];
	double y;
#endif
	double r;

	mx2 = mx[2];
	my2 = my[2];

	/* exceptions: x or y = {0,inf,NaN} */
	if (mx2 == 0 || my2 == 0 || mx2 >= 0x7ff00000 || my2 >= 0x7ff00000) {
		mz[0] = (mx[0] ^ my[0]) | 1;
		mz[1] = 0;
		for (i = 3; i < nw; i++)
			mz[i] = 0;
		if (mx2 == 0) {
			if (my2 > 0x7ff00000)
				mz[2] = my2;
			else /* 0/NaN is NaN */
			if (my2 == 0)
				mz[2] = 0x7fffffff;
			else /* 0/0 is sNaN */
				mz[2] = mx2; /* 0/x is 0 */
		} else if (mx2 == 0x7ff00000) {	/* inf/NaN is NaN */
			if (my2 > mx2)
				mz[2] = my2;
			else /* inf/inf is sNaN */
			if (my2 == mx2)
				mz[2] = 0x7fffffff;
			else	/* inf/x is inf */
				mz[2] = mx2;
		} else {	/* mx2 > 0x7ff00000 */
			if (my2 > mx2)
				mz[2] = my2;	/* NaN/NaN is NaN */
			else
				mz[2] = mx2;
		}
		return;
	}
	nwp = nw + 2;		/* extra words precision */

	myi = (int *) calloc(nwp, sizeof (int));
	mri = (int *) calloc(nwp, sizeof (int));
	mq = (int *) calloc(nwp, sizeof (int));
	mi_two = (int *) calloc(nwp, sizeof (int));

	/* initial value of mri ~ 1/|my| */
	mri[0] = 1;
	if (((my[2] - 1) | my[3] | my[4]) == 0) {
		mri[1] = -my[1];
		mri[2] = 1;
		for (i = 3; i < nwp; i++)
			mri[i] = 0;
	} else {
		mri[1] = -1 - my[1];
		if (my[4] != 0)
			r = twom24 * (double) my[4];
		else
			r = 0.0;
		if (my[3] != 0)
			r += (double) my[3];
		r = two24 / ((double) my[2] + twom24 * r);
		mri[2] = (int) r;
		r = (r - (double) mri[2]) * two24;
		mri[3] = (int) r;
		mri[4] = (int) ((r - (double) mri[3]) * two24);
		for (i = 5; i < nwp; i++)
			mri[i] = 0;
	}
#ifdef DEBUG
	y = __libmcr_mi_mitod(my, nw, mtail);
	printf("y   = % 1.20e %08X %08X\n", y, y);
	r = __libmcr_mi_mitod(mri, nw, mtail);
	printf("r   = % 1.20e %08X %08X\n", r, r);
	r = 1.0 / y;
	printf("1/y = % 1.20e %08X %08X\n", r, r);
#endif

	/*
	 * myi = my mq  = mx loop until r=1+tiny, where tiny is less than 1/2
	 * of the mantissa y = y*r q = q*r r = 2-y	r is slightly over
	 * 1.0 Finally return q*r (the last update);
	 */
	mi_two[0] = 1;
	mi_two[1] = 0;
	mi_two[2] = 2;
	for (j = 3; j < nwp; j++)
		mi_two[j] = 0;
	for (j = 1; j < nw; j++) {
		myi[j] = my[j];
		mq[j] = mx[j];
	}
	for (j = nw; j < nwp; j++) {
		myi[j] = 0;
		mq[j] = 0;
	}
	myi[0] = mq[0] = 1;

	i = 1;
	while (i) {
		__libmcr_mm_mul(myi, mri, myi, nwp);
		__libmcr_mm_mul(mq, mri, mq, nwp);
		__libmcr_mm_sub(mi_two, myi, mri, nwp);
#ifdef DEBUG
		k = 2;
		if (mri[1] == 0 && mri[2] == 1)
			for (k = 3; mri[k] == 0 && k < nw; k++);
		k -= 2;
		printf("loop: r has %d zero after 1. (%d percent zeros)\n",
			k, (int) ((100.0 * k) / (double) nw));
#endif
		if (mri[1] == 0) {
			i = mri[2] - 1;
			k = ((nwp) >> 1) + 1;
			for (j = 3; j <= k; j++)
				i |= mri[j];
			if (i != 0)
				i = 1;
		}
	}
	__libmcr_mm_mul(mq, mri, mq, nwp);

	/* restore sign */
	mq[0] = (mx[0] ^ my[0]) | 1;

	for (j = 0; j < nw; j++) {
		mz[j] = mq[j];
	}
#ifdef DEBUG
	for (j = 0; j <= nw; j++) {
		printf("  mq[%3d] = %08X\n", j, mq[j]);
	}
#endif
	free(mq);
	free(myi);
	free(mri);
	free(mi_two);
}
