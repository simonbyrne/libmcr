
#pragma ident "@(#)__libmcr_k_mi_rem_pio2.c 1.4 04/02/25 SMI"

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
 * int k_mi_rem_pio2(double x, int* mc, int nw)
 *
 * This routine accurately compute the mathematical value of
 *	mc = x - N*(2/pi),
 * where N is define to be the (mathematical) integer nearest
 * to x*(2/pi), and return n = N mod 8.
 *
 * Method:
 *	See "Argument Reduction for Huge Arguments: Good to the Last Bit"
 *	by K.C. Ng, June 1992 (unpublished).
 *
 * Implementation notes:
 *    1 Let two_pi[] be the mi number of 2/pi. Then each integer part of
 *	two_pi[k] is corresponding to the actual value
 *		two_pi[k] * 2**(24*(1-k))	 ... assuming k>=2
 *    2 Let m = ilogb(x) (the unbias exponent of x). Then the value of x
 *	can be written as
 *		2**(m-52) * X
 *	where X is significant part of x that is regard as an integer.
 *	Thus, if
 *		2**(m-52) * X  *  two_pi[k] * 2**(24*(1-k))
 *		-------------------------------------------
 *				    8
 *
 *    		= 2**[(m-52)-3+24*(1-k)]*(X*two_pi[k])
 *      (*)	= 2**[m-31 - 24*k]  * (X*two_pi[k])
 *	is an integer, then the first k+1 word of two_pi can be regarded as
 *	the leading part of two_pi such that the multiplication of it with x
 *	is equal to zero mod 8.
 *
 *	For (*) to be an integer is equivalent to find the k such that
 *		m-31-24*k is >= 0   <==>   (m-31)/24 >= k
 *	Thus, let k  =  (m-31)/24 (truncated to an integer). If k<2, then
 *	no need to replace two_pi. Otherwise, let
 *		two_pi_trailing[0] = two_pi[0]
 *		two_pi_trailing[1] = -k
 *		two_pi_trailing[j] = two_pi[k+j-1]    ...for j>=2
 *	The we can replace two_pi by two_pi_trailing in the computation of
 *	the argument reduction.
 *
 *    3 Since the worst case of double produce > 2**-66 reduced argument,
 *	so three more words will suffice to compensate for the tiny
 *	reduction. Allow one more word for the loss in choosing k in step 2,
 *	then we may use nw+4 to guarantee full accuracy of argument reduction.
 */

#ifdef DEBUG
#define	H6(x)	x[0], x[1], x[2], x[3], x[4], x[5]
#define	HA(x)	*(int *)&x, *(1+(int *)&x)
#include <stdio.h>
#endif
#include <stdlib.h>
#include "mcr.h"

extern const long _TBL_ipio2[];	/* constant up to 66 words */
extern const long _TBL_pio2[];	/* constant up to nw=52 */

int
__libmcr_k_mi_rem_pio2(x, mc, nw)
	double x;
	int mc[], nw;
{
	int *mx, *mr, *mt, *two_pi, *pio2, i, j, k, hx, ix, nwp, N, nnw,
	*ipio2, *mi_one;

	hx = HIGH_WORD(x);
	ix = hx & 0x7fffffff;
	/* return x if |x|<pi/4 */
	if (ix <= 0x3fe921fb) {
		if (ix < 0x3fe921fb ||
			(((unsigned) LOW_WORD(x)) <= (unsigned) 0x54442d18)) {
			__libmcr_mi_dtomi(x, mc, nw);
			return (0);
		}
	}
	/* now |x|>pi/4 */

	nwp = nw + 4;
	mx = (int *) calloc(nwp, sizeof (int));
	mt = (int *) calloc(nwp, sizeof (int));
	mr = (int *) calloc(nwp, sizeof (int));
	two_pi = (int *) calloc(nwp, sizeof (int));
	pio2 = (int *) calloc(nwp, sizeof (int));
	__libmcr_mm_pio2(pio2, nwp);

	/* setup two_pi = 2/pi in mi format */
	k = ((ix >> 20) - 1023 - 31) / 24;
#ifdef DEBUG
	printf("ix     = %08X\n", ix);
	printf("ix>>20 = %08X\n", ix >> 20);
	printf("k      = %d\n", k);
#endif
	if (k < 2)
		k = 1;
	two_pi[0] = 1;
	two_pi[1] = -k;
	if (nwp + k - 3 < 66) {
		for (i = 2; i < nwp; i++)
			two_pi[i] = __libmcr_TBL_ipio2[i - 3 + k];
#ifdef DEBUG
		printf("nwp  = %d, using constant\n", nwp);
		for (i = 0; i < nwp; i++)
			printf("%08X\n", two_pi[i]);
#endif
	} else {
		nnw = nwp + k;
		ipio2 = (int *) calloc(nnw, sizeof (int));
		mi_one = (int *) calloc(nnw, sizeof (int));
		for (i = 1; i < nnw; i++)
			mi_one[i] = 0;
		mi_one[0] = mi_one[2] = 1;
		__libmcr_mm_pio2(ipio2, nnw);
		__libmcr_mm_div(mi_one, ipio2, ipio2, nnw);
		for (i = 2; i < nwp; i++)
			two_pi[i] = ipio2[i - 1 + k];
		free(ipio2);
		free(mi_one);
#ifdef DEBUG
		printf("nnw  = %d, using __libmcr_mm_pio2 alogrithm\n", nnw);
		for (i = 0; i < nwp; i++)
			printf("%08X\n", two_pi[i]);
#endif
	}

#ifdef DEBUG
	printf("x = %08X %08X\n", HA(x));
#endif
	/* |x| * 2/pi */
	__libmcr_mi_dtomi(x, mx, nwp);
	mx[0] = 1;
	__libmcr_mm_mul(mx, two_pi, mr, nwp);

#ifdef DEBUG
	printf("x* short 2/pi =\n");
	for (i = 0; i < nwp; i++)
		printf("%08X\n", mr[i]);
#endif
	/* ignore part that is >= 2**24 */
	i = mr[1];
	if (i > 0) {
		mr[1] = 0;
		for (j = 2; j < nwp - i; j++)
			mr[j] = mr[j + i];
		for (j = nwp - i; j < nwp; j++)
			mr[j] = 0;
	}
#ifdef DEBUG
	printf("after adjustment:\n");
	for (i = 0; i < nwp; i++)
		printf("%08X\n", mr[i]);
#endif

	/* determine N = x*2/pi rounded to an integer */
	N = 0;
	if (i == -1) {
		if (mr[2] >= 0x800000)
			N = 1;
	} else {
		N = mr[2];
		if (mr[3] >= 0x800000)
			N += 1;
	}

#ifdef DEBUG
	printf("N = %d %08X\n", N, N);
#endif

	/* compute mc = [ (x*2/pi) - N ] * pi/2 */
	mt[0] = -1;
	mt[1] = 0;
	mt[2] = N;
	for (i = 3; i < nwp; i++)
		mt[i] = 0;
#ifdef DEBUG
	printf("N in mi = %08X %08X %08X %08X %08X %08X\n", H6(mt));
#endif
	__libmcr_mm_add(mr, mt, mt, nwp);
#ifdef DEBUG
	printf("(x*2/pi)-N  = %08X %08X %08X %08X %08X %08X\n", H6(mt));
#endif
	__libmcr_mm_mul(mt, pio2, mc, nw);
#ifdef DEBUG
	printf("pi/2*((x*2/pi)-N)  = %08X %08X %08X %08X %08X %08X\n", H6(mc));
#endif
	if (hx < 0) {
		mc[0] = -mc[0];
		N = -N;
	}
	free(pio2);
	free(mx);
	free(mt);
	free(mr);
	free(two_pi);
	return (N & 7);
}
