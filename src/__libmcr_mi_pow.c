
#pragma ident "@(#)__libmcr_mi_pow.c 1.6 04/02/25 SMI"

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
 * double __libmcr_mi_pow(double x, double y, int nw, int *ncrd, int rnd)
 *
 * Using multi-precision arithmetic, __libmcr_mi_pow return an almost
 * correctly rounded pow(x,y), with *ncrd=1 indicate if it it not guarantee
 * correctly rounded. Here nw is the number of integer words in the
 * multi-precision array. (Three words used for sign, exponent, and
 * leading word. Each addition integer represented an additional
 * 24-bits precision.)
 *
 * Method:
 *	pow(x,y) = exp(y*log(x)).
 *	The argument reduction step for exp on y*log(x) may involve
 *	cancellation:  r = y*log(x) - n*ln2. We should start with
 *	four extra words to guard the cancellation, and check. If
 *	not enough, we will increase the accuracy.
 *
 *
 * Special cases:
 *	1.  (anything) ** 0  is 1
 *	2.  (anything) ** 1  is itself
 *	3.  (anything) ** NAN is NAN
 *	4.  NAN ** (anything except 0) is NAN
 *	5.  +-(|x| > 1) **  +INF is +INF
 *	6.  +-(|x| > 1) **  -INF is +0
 *	7.  +-(|x| < 1) **  +INF is +0
 *	8.  +-(|x| < 1) **  -INF is +INF
 *	9.  +-1         ** +-INF is NAN
 *	10. +0 ** (+anything except 0, NAN)               is +0
 *	11. -0 ** (+anything except 0, NAN, odd integer)  is +0
 *	12. +0 ** (-anything except 0, NAN)               is +INF
 *	13. -0 ** (-anything except 0, NAN, odd integer)  is +INF
 *	14. -0 ** (odd integer) = -( +0 ** (odd integer) )
 *	15. +INF ** (+anything except 0,NAN) is +INF
 *	16. +INF ** (-anything except 0,NAN) is +0
 *	17. -INF ** (anything)  = -0 ** (-anything)
 *	18. (-anything) ** (integer) is (-1)**(integer)*(+anything**integer)
 *	19. (-anything except 0 and inf) ** (non-integer) is NAN
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
#ifdef DEBUG
#include <stdio.h>
#define	H2(x) *(int *)&x, *(1+(int *)&x)
#endif

static const double
half[] = {0.5, -0.5},
/* Redundant exception handling code: see __libmcr_pow */
#if 0
huge = 1.0e300,
zero = 0.0,
one = 1.0,
tiny = 1.0e-300,
#endif
invln2 = 1.44269504088896338700e+00,	/* 0x3ff71547, 0x652b82fe */
two54 = 1.80143985094819840000e+16;

double
__libmcr_mi_pow(x, y, nw, ncrd, rnd)
	double x, y;
	int nw, *ncrd, rnd;
{
	/*
	 * rnd: rounding direction: 0 --  nearest, 1 - chopped, 2 - up, 3 -
	 * down
	 */
	int k, n, *mc, *mt, mtail[4], *my, *ln2;
	int hx, ix, na, nwp;
	double z, ax, w;
/* Redundant exception handling code: see __libmcr_pow */
#if 0
	int iy, hy, ly, lx, j, yisint;
#endif


	/* determine exception cases */
	*ncrd = 0;

	hx = HIGH_WORD(x);
	ix = hx & 0x7fffffff;
#ifdef DEBUG
	printf(" x = %08X %08X %1.20e\n", H2(x), x);
#endif

/* Redundant exception handling code: see __libmcr_pow */
#if 0
	hy = HIGH_WORD(y);
	lx = LOW_WORD(x);
	ly = LOW_WORD(y);
	iy = hy & 0x7fffffff;
	/* y==zero: x**0 = 1 */
	if ((iy | ly) == 0)
		return (one);

	/* +-NaN return x+y */
	if (ix > 0x7ff00000 || ((ix == 0x7ff00000) && (lx != 0)) ||
	    iy > 0x7ff00000 || ((iy == 0x7ff00000) && (ly != 0)))
		return (x + y);

	/*
	 * determine if y is an odd int when x < 0 yisint = 0	... y is not
	 * an integer yisint = 1	... y is an odd int yisint = 2	... y
	 * is an even int
	 */
	yisint = 0;
	if (hx < 0) {
		if (iy >= 0x43400000)
			yisint = 2;	/* even integer y */
		else if (iy >= 0x3ff00000) {
			k = (iy >> 20) - 0x3ff;	/* exponent */
			if (k > 20) {
				j = ly >> (52 - k);
				if ((j << (52 - k)) == ly)
					yisint = 2 - (j & 1);
			} else if (ly == 0) {
				j = iy >> (20 - k);
				if ((j << (20 - k)) == iy)
					yisint = 2 - (j & 1);
			}
		}
	}
	/* special value of y */
	if (ly == 0) {
		if (iy == 0x7ff00000) {	/* y is +-inf */
			if (((ix - 0x3ff00000) | lx) == 0)
				return (y - y);	/* inf**+-1 is NaN */
			else if (ix >= 0x3ff00000)
			{
				/* (|x|>1)**+-inf = * inf,0 */
				return ((hy >= 0) ? y : zero);
			} else	/* (|x|<1)**-,+inf = inf,0 */
				return ((hy < 0) ? -y : zero);
		}
		if (iy == 0x3ff00000) {	/* y is  +-1 */
			if (hy < 0)
				return (one / x);
			else
				return (x);
		}
		if (hy == 0x40000000)
			return (x * x);	/* y is  2 */
		if (hy == 0x3fe00000) {	/* y is  0.5 */
			if (hx >= 0)	/* x >= +0 */
				return (sqrt(x));
		}
	}
	ax = fabs(x);
	/* special value of x */
	if (lx == 0) {
		if (ix == 0x7ff00000 || ix == 0 || ix == 0x3ff00000) {
			z = ax;	/* x is +-0,+-inf,+-1 */
			if (hy < 0)
				z = one / z;	/* z = (1/|x|) */
			if (hx < 0) {
				if (((ix - 0x3ff00000) | yisint) == 0) {
					/* (-1)**non-int is NaN */
					z = (z - z) / (z - z);
				} else if (yisint == 1)
					z = -z;	/* (x<0)**odd = -(|x|**odd) */
			}
			return (z);
		}
	}
	/* (x<0)**(non-int) is NaN */
	if ((((hx >> 31) + 1) | yisint) == 0)
		return ((x - x) / (x - x));

	/* |y| is huge */
	if (iy > 0x41e00000) {	/* if |y| > 2**31 */
		if (iy > 0x43f00000) {	/* if |y| > 2**64, must o/uflow */
			if (ix <= 0x3fefffff)
				return ((hy < 0) ? huge * huge : tiny * tiny);
			if (ix >= 0x3ff00000)
				return ((hy > 0) ? huge * huge : tiny * tiny);
		}
		/* over/underflow if x is not close to one */
		if (ix < 0x3fefffff)
			return ((hy < 0) ? huge * huge : tiny * tiny);
		if (ix > 0x3ff00000)
			return ((hy > 0) ? huge * huge : tiny * tiny);
	}
#else
	ax = fabs(x);
#endif

	/* compute log(x) */
	nwp = nw + 4;
	mc = (int *) calloc(nwp, sizeof (int));
	mt = (int *) calloc(nwp, sizeof (int));
	my = (int *) calloc(nwp, sizeof (int));
	ln2 = (int *) calloc(nwp, sizeof (int));
	__libmcr_mm_ln2(ln2, nwp);

	na = 0;

	if (ix < 0x00100000) {	/* subnormal */
		x *= two54;
		na = -54;
		ix = HIGH_WORD(x);
	}
	k = ((ix + 0x95f61) >> 20) - 0x3ff;
	na += k;
	ix -= (k << 20);	/* normalize x to [1/sqrt2, sqrt2] */
	HIGH_WORD(x) = ix;

#ifdef DEBUG
	printf(" k = %d,  x = %08X %08X %1.20e\n", k, H2(x), x);
#endif

	/* Log(x) = n*ln2 + log(r) */
	__libmcr_k_mi_log(x, mc, nwp);
	__libmcr_mm_muli(ln2, na, mt, nwp);
	__libmcr_mm_add(mt, mc, mc, nwp);

	/* w = y*log(x) */
	__libmcr_mi_dtomi(y, my, nwp);
	__libmcr_mm_mul(my, mc, mc, nwp);
	w = __libmcr_mi_mitod(mc, nwp, mtail);

	/* n =  w*invln2 + sign(0.5,w) */
	n = w * invln2 + half[(((unsigned) mc[0]) >> 31)];

	/* loop if cancellation occurs on y*log(x) - n*ln2 */
	__libmcr_mm_muli(ln2, n, my, nwp);
	__libmcr_mm_sub(mc, my, mc, nwp);
	k = 93;
#ifdef DEBUG
	printf("my = %08X %08X %08X %08X\n", my[0], my[1], my[2], my[3]);
	printf("mc = %08X %08X %08X %08X\n", mc[0], mc[1], mc[2], mc[3]);
	printf("mi_ilogb(my) : %08X\n", __libmcr_mi_ilogb(my));
	printf("mi_ilogb(mc) : %08X\n", __libmcr_mi_ilogb(mc));
	fflush(stdout);
#endif
	while (__libmcr_mi_ilogb(my) - __libmcr_mi_ilogb(mc) >= k) {
#ifdef DEBUG
		printf("nwp,k = %d %d\n", nwp, k);
		printf("my= %08X %08X %08X %08X\n", my[0], my[1], my[2], my[3]);
		printf("mc= %08X %08X %08X %08X\n", mc[0], mc[1], mc[2], mc[3]);
		printf("mi_ilogb(my) : %08X\n", __libmcr_mi_ilogb(my));
		printf("mi_ilogb(mc) : %08X\n", __libmcr_mi_ilogb(mc));
		fflush(stdout);
#endif
		nwp += 2;
		k += 48;
		free(ln2);
		free(mc);
		free(mt);
		free(my);
		ln2 = (int *) calloc(nwp, sizeof (int));
		mc = (int *) calloc(nwp, sizeof (int));
		mt = (int *) calloc(nwp, sizeof (int));
		my = (int *) calloc(nwp, sizeof (int));
		__libmcr_mm_ln2(ln2, nwp);
		__libmcr_k_mi_log(x, mc, nwp);
		__libmcr_mm_muli(ln2, na, mt, nwp);
		__libmcr_mm_add(mt, mc, mc, nwp);
		__libmcr_mi_dtomi(y, my, nwp);
		__libmcr_mm_mul(my, mc, mc, nwp);
		w = __libmcr_mi_mitod(mc, nwp, mtail);
		n = w * invln2 + half[(((unsigned) mc[0]) >> 31)];
		__libmcr_mm_muli(ln2, n, my, nwp);
		__libmcr_mm_sub(mc, my, mc, nwp);
	}

	/* exp(r) */
	__libmcr_k_mm_exp(mc, mc, nw);

	/* scalb mc by 2**n */
	__libmcr_mm_scalbn(mc, nw, n);

	/* final rounding */
	z = __libmcr_mi_final(mc, nw, ncrd, rnd);
	free(ln2);
	free(mc);
	free(mt);
	free(my);
	return (z);
}
