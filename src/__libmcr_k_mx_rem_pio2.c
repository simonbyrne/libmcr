
#pragma ident "@(#)__libmcr_k_mx_rem_pio2.c 1.6 04/02/25 SMI"

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
 * int __libmcr_k_mx_rem_pio2(double x, double y[])
 *
 * return the remainder of x rem pi/2 in y[0]+y[1]
 * use __libmcr_k_rem_pio2m() for really huge x
 */

#ifdef debug
#include <stdio.h>
#define	PNT(x) printf("% 1d  % 4d  %x.%06X%06X%06X%06X%06X\n", \
	x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7])
#endif
#include "mcr.h"

/*
 * INVPIO2:  53 bits of 2/pi
 * PIO2_1:   first  33 bit of pi/2
 * PIO2_1T:  pi/2 - PIO2_1
 * PIO2_2[0]:second 33 bit of pi/2
 * PIO2_2[1]:pi/2 - PIO2_1 - PIO2_2[0]
 * PIO2_2[2]:third 33 bit of pi/2
 * PIO2_2[3]:pi/2 - PIO2_1 - PIO2_2[0] - PIO2_2[2]
 * PIO2_2[4]:fourth 33 bit of pi/2
 * PIO2_2[5]:pi/2 - PIO2_1 - PIO2_2[0] - PIO2_2[2] - PIO2_2[4]
 * PIO2_2[6]:fifth 33 bit of pi/2
 * PIO2_2[7]:pi/2 - PIO2_1 - PIO2_2[0] - PIO2_2[2] - PIO2_2[4] - PIO2_2[6]
 * PIO2_2[8]:sixth 20 bit of pi/2
 * PIO2_2[9]:pi/2 -PIO2_1 -PIO2_2[0] -PIO2_2[2] -PIO2_2[4] -PIO2_2[6] -PIO2_2[8]
 */
static const double rm[] = {
	/* PIO2_2  = */
	6.07710050630396597660e-11,	/* in hex: 3DD0B461, 1A600000 */
	2.02226624879595063154e-21,	/* in hex: 3BA3198A, 2E037073 */
	2.02226624871116645580e-21,	/* in hex: 3BA3198A, 2E000000 */
	8.47842766036889956997e-32,	/* in hex: 397B839A, 252049C1 */
	8.47842766034822932286e-32,	/* in hex: 397B839A, 25200000 */
	2.06703210982639879819e-43,	/* in hex: 37127044, 533E63A0 */
	2.06703210945082419395e-43,	/* in hex: 37127044, 53300000 */
	3.75574629696355818273e-53,	/* in hex: 350CC740, 20BBEA63 */
	3.75574629660149376100e-53,	/* in hex: 350CC740, 20B00000 */
	3.62064742712430339369e-63,	/* in hex: 32F7D4C7, 6273644A */
	/* INVPIO2 = */
	0.636619772367581343075535,	/* 2^ -1  * Hex * 1.45F306DC9C883 */
	/* PIO2_1  = */
	1.57079632673412561417,
	/* PIO2_1T = */
	6.07710050650619260148e-11,
	/* HALF    = */
	0.5,
};

#define	PIO2_2	rm
#define	INVPIO2	rm[10]
#define	PIO2_1	rm[11]
#define	PIO2_1T	rm[12]
#define	HALF	rm[13]

int
__libmcr_k_mx_rem_pio2(double x, double y[])
{
	double z, w, t, r, fn, y0;
	double tx[3];
	int e0, i, j, nx, n, ix, hx, lx, hr, i0;
#ifdef debug
	int ms[10], mw[10], mx[10], my[10], mz[10], nw = 8, mc[10], nz,
		mt[10], mr[10], my0[10];

	nz = __libmcr_k_mi_rem_pio2(x, my, nw);
#endif

	z = fabs(x);
	w = z * INVPIO2;
	hx = HIGH_WORD(x);
	lx = LOW_WORD(x);
	ix = hx & 0x7fffffff;
	/* for |x| <= 1638400 */
	if (ix <= 0x41390000) {
		n = (int) (w + HALF);	/* |x| ~ n*pi/2, n < 2^19 */
#ifdef debug
	printf(" w = %1.20e %08X %08X \n", w, .HIGH_WORD(w), .LOW_WORD(w));
	printf(" HALF = %1.20e\n", HALF);
	printf(" n = %d\n", n);
#endif
		if (ix <= 0x3fe921fa) {
			y[0] = x;
			y[1] = 0;
			return (0);
		}		/* |x| < pi/2 */
		fn = (double) n;
		j = ix >> 20;
		t = fn * PIO2_1;
		w = fn * PIO2_1T;	/* 1st round good to 85 bit */
		r = z - t;
		i0 = 0;
		y0 = r - w;
		hr = HIGH_WORD(r);
		i = (hr & (~0x80000000)) >> 20;
#ifdef debug
		printf("x        = %08X %08X % 1.20e\n",
			HIGH_WORD(x), LOW_WORD(x), x);
		printf("n        = %d\n", n);
		printf("nz       = %d\n", nz);
		printf("n*PIO2_1 = %08X %08X % 1.20e\n",
			HIGH_WORD(t), LOW_WORD(t), t);
		printf("n*PIO2_1T= %08X %08X % 1.20e\n",
			HIGH_WORD(w), LOW_WORD(w), w);
		__libmcr_mi_dtomi(t, mt, nw);
		__libmcr_mi_dtomi(w, mw, nw);
		__libmcr_mi_dtomi(z, mz, nw);
		__libmcr_mm_sub(mz, mt, mr, nw);
		__libmcr_mm_sub(mr, mw, my0, nw);
		printf("z - n*PIO2 =\n");
		PNT(my);
		PNT(my0);
		printf("r = z-t  = %08X %08X % 1.20e\n",
			HIGH_WORD(r), LOW_WORD(r), r);
		printf("y0= r-w  = %08X %08X % 1.20e\n",
			HIGH_WORD(y0), LOW_WORD(y0), y0);
		printf("j = %d, i = %d, j-i=%d > 6?\n", j, i, j - i);
#endif
		while ((j - i) >= 7) { /* 7 bit cancellation (insure 78 bits) */
			t = fn * PIO2_2[i0];  /* at most five loop is suffice */
			z = r;
			j = i;
#ifdef debug
			printf("Are r-t exact?\n");
			printf("r = %08X %08X %1.20e\n", r, r);
			printf("t = %08X %08X %1.20e\n", t, t);
			__libmcr_mi_dtomi(t, mt, nw);
			__libmcr_mi_dtomi(r, mw, nw);
			__libmcr_mm_sub(mw, mt, mt, nw);
#endif
			r -= t;
#ifdef debug
			__libmcr_mi_dtomi(r, mw, nw);
			__libmcr_mm_sub(mw, mt, mt, nw);
			if (mt[2] == 0)
				printf("Yes.\n");
			else
				printf("No.\n");
#endif
			w = fn * PIO2_2[i0 + 1];
			y0 = r - w;
			i0 += 2;
			hr = HIGH_WORD(r);
			i = (hr & (~0x80000000)) >> 20;
#ifdef debug
			printf("w = %08X %08X %1.20e\n", w, w);
			printf("y0= %08X %08X %1.20e\n", y0, y0);
			__libmcr_mi_dtomi(t, mt, nw);
			__libmcr_mi_dtomi(w, mw, nw);
			__libmcr_mm_sub(mr, mt, mr, nw);
			__libmcr_mm_sub(mr, mw, my0, nw);
			printf("z - n*PIO2[%1d] =\n", i0 + 2);
			PNT(my);
			PNT(my0);
			printf("new j = %d, i = %d, j-i=%d > 5?\n",
				j, i, j - i);
#endif
		}
		if (i0 == 0 || (j - i) > 1) {
			r -= y0;
			y[0] = y0;
			y[1] = r - w;
		} else {
			t += y0 - z;
			y[0] = y0;
			y[1] = -w - t;
		}

#ifdef debug
		printf("test final sum (n=%d):\n", n);
		__libmcr_mi_dtomi(y[0], mr, nw);
		printf("y[0]  = % 1.20e ", y[0]);
		PNT(mr);
		__libmcr_mi_dtomi(y[1], mt, nw);
		printf("y[1]  = % 1.20e ", y[1]);
		PNT(mt);
		__libmcr_mm_add(mr, mt, mt, nw);
		printf("y0+y1 = % 1.20e ", y[0] + y[1]);
		PNT(mt);
		PNT(mt);
#endif

		if (hx < 0) {
			y[1] = -y[1];
			y[0] = -y0;
			return (-n);
		} else {
			return (n);
		}
	}
	if (ix >= 0x7ff00000) {	/* x is inf or NaN */
		y[0] = y[1] = x - x;
		return (0);
	}
	e0 = (ix >> 20) - 1046;	/* e0 = ilogb(z)-23; */
	/* break z into three 24 bit pieces */
	i = ((lx & 0x1f) << 19);
	tx[2] = (double) i;
	j = ((lx >> 5) & 0xffffff);
	tx[1] = (double) j;
	tx[0] = (double) ((((ix & 0xfffff) | 0x100000) << 3) |
			((unsigned) lx >> 29));
	nx = 3;
	if (i == 0) {
		nx--;
		if (j == 0)
			nx--;
	}			/* skip zero term */
	n = __libmcr_k_rem_pio2m(tx, y, e0, nx, 2);
	if (hx < 0) {
		y[0] = -y[0];
		y[1] = -y[1];
		return (-n);
	}
	return (n);
}
