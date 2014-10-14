
#pragma ident "@(#)__libmcr_mm_sqrt.c 1.3 04/02/25 SMI"

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
 * mm_sqrt(mx,mr,nw)
 * Compute mr = sqrt(mx) in mi precision, nw the number of words.
 *
 * mz = ~1/sqrt(mx)
 * 	mr = mx*mz
 * 	my = mx*mz*mz
 * 	mz = (3-my)/2
 * Loop mz = 1+tiny, where tiny is less than 1/2 of the mantissa
 *	my = my*mz*mz
 *	mr = mr*mz
 *	mz = (3-my)/2
 * Finally,
 *	return mr = mr*mz.
 *
 * Analysis:
 *	z0 ~ 1/sqrt(x), r0 = x*z0 ~ sqrt(x),  y0 = x*z0*z0 ~ 1
 *	z1 = (3-y0)/2
 * For i = 0, 1, 2,...
 *	y(i+1) = yi*z(i+1)*z(i+1)
 *	r(i+1) = ri*z(i+1)
 *	z(i+2) = (3-y(i+1))/2
 *
 * If z0 = 1/sqrt(x)*(1-E0), then
 *	y0 = (1-E0)^2;
 *	r0 = sqrt(x)*(1-E0);
 *	z1 = (3-(1-E0)^2)/2 = (2+2E0-E0^2)/2 = 1+ E0-E0^2/2
 *
 *	y1 = y0*z1*z1 = (1-E0)^2*(1+E0-E0^2/2)^2 = [(1-E0)(1+E0-E0^2/2)]^2
 *	   = [1+E0-E1-E0^2/2-E0^2+E0^3/2]^2 = [1-E0^2/2+E0^3/2]^2
 *	   == [1-E1]^2, where E1 = E0^2/2-E0^3/2
 *	r1 = r0*z1 = sqrt(x)*(1-E0)*(1+E0-E0^2/2) = sqrt(x)*(1-E1)
 *	z2 = (3-y1)/2 = (3-(1-E1)^2)/2 = 1 + E1 - E1^2/2
 *
 *	y2 = [1-E2]^2, where E2 = E1^2/2 - E1^3/2
 *	r2 = r1*z2 = sqrt(x)*(1-E1)*(1+E1-E1^2/2) = sqrt(x)*(1-E2)
 *	z2 = (3-y2)/2 = (3-(1-E2)^2)/2 = 1 + E2 - E2^2/2
 *
 *	...
 *
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
__libmcr_mm_sqrt(mx, mr, nw)
	int mx[], mr[], nw;
{
	int i, k, j, *myi, *mzi, *mri, *mi_three, nwp;
	double y, r;
#ifdef DEBUG
	int mtail[4];
	double x;
#endif

	/* exception sqrt(+-0)=+-0, sqrt(-ve) = NaN */
	if (mx[2] == 0) {
		mr[0] = mx[0];
		for (i = 1; i < nw; i++)
			mr[i] = 0;
		return;
	}
	if (mx[0] < 0) {
		mr[0] = mx[0];
		mr[1] = 0;
		mr[2] = 0x7fffffff;
		for (i = 3; i < nw; i++)
			mr[i] = 0;
		return;
	}
	nwp = nw + 2;
	myi = (int *) calloc(nwp + 1, sizeof (int));
	mri = (int *) calloc(nwp + 1, sizeof (int));
	mzi = (int *) calloc(nwp + 1, sizeof (int));
	mi_three = (int *) calloc(nwp + 1, sizeof (int));

	/* initial valure mzi ~ 1/sqrt(x) */

	mzi[0] = 1;
	mzi[1] = -(mx[1] >> 1) - 1;
	if ((mx[1] & 1) == 0)
		y = two24;
	else
		y = 4096.0;

	if (mx[4] != 0)
		r = twom24 * (double) mx[4];
	else
		r = 0.0;
	if (mx[3] != 0)
		r += (double) mx[3];
	r = y / sqrt((double) mx[2] + twom24 * r);
	mzi[2] = (int) r;
	r = (r - (double) mzi[2]) * two24;
	mzi[3] = (int) r;
	mzi[4] = (int) ((r - (double) mzi[3]) * two24);
	for (i = 5; i < nwp; i++)
		mzi[i] = 0;


#ifdef DEBUG
	x = __libmcr_mi_mitod(mx, nw, mtail);
	printf("x         = % 1.20e %08X %08X\n", x, x);
	r = __libmcr_mi_mitod(mzi, nw, mtail);
	printf("r         = % 1.20e %08X %08X\n", r, r);
	r = 1.0 / sqrt(x);
	printf("1/sqrt(x) = % 1.20e %08X %08X\n", r, r);
#endif

	/*
	 * mz = ~1/sqrt(mx) mr = mx*mz my = mx*mz*mz	(or my = mr*mz) mz =
	 * (3-my)/2 Loop mz = 1+tiny, where tiny is less than 1/2 of the
	 * mantissa my = my*mz*mz mr = mr*mz mz = (3-my)/2 Finally, return mr
	 * = mr*mz.
	 */
	mi_three[0] = 1;
	mi_three[1] = 0;
	mi_three[2] = 3;
	for (j = 3; j < nwp; j++)
		mi_three[j] = 0;
	for (j = 1; j < nw; j++)
		myi[j] = mx[j];
	for (j = nw; j < nwp; j++)
		myi[j] = 0;
	myi[0] = 1;

	__libmcr_mm_mul(myi, mzi, mri, nwp);
	__libmcr_mm_mul(mzi, mri, myi, nwp);
	__libmcr_mm_sub(mi_three, myi, mzi, nwp);
	__libmcr_mm_scalbn(mzi, nwp, -1);

	i = 1;
	while (i) {
		__libmcr_mm_mul(myi, mzi, myi, nwp);
		__libmcr_mm_mul(myi, mzi, myi, nwp);
		__libmcr_mm_mul(mri, mzi, mri, nwp);
		__libmcr_mm_sub(mi_three, myi, mzi, nwp);
		__libmcr_mm_scalbn(mzi, nwp, -1);
#ifdef DEBUG
		k = 2;
		if (mzi[1] == 0 && mzi[2] == 1)
			for (k = 3; mzi[k] == 0 && k < nw; k++);
		k -= 2;
		printf("loop: r has %d zero after 1. (%d percent zeros)\n",
			k, (int) ((100.0 * k) / (double) nw));
#endif
		if (mzi[1] == 0) {
			i = mzi[2] - 1;
			k = ((nwp) >> 1) + 1;
			for (j = 3; j <= k; j++)
				i |= mzi[j];
			if (i != 0)
				i = 1;
		}
	}
	__libmcr_mm_mul(mri, mzi, mri, nwp);

	for (j = 0; j < nw; j++) {
		mr[j] = mri[j];
	}
#ifdef DEBUG
	for (j = 0; j <= nw; j++) {
		printf("  mri[%3d] = %08X\n", j, mri[j]);
	}
#endif
	free(mzi);
	free(myi);
	free(mri);
	free(mi_three);
}
