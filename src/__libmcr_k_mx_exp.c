
#pragma ident "@(#)__libmcr_k_mx_exp.c 1.6 04/02/25 SMI"

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
 * double __libmcr_k_mx_exp((double)x,(double*)e,(int*)m)
 *
 * Multi-precision exp with correction term e up to 22 extra bits.
 * To avoid over/underflow issue, the result of exp(x) is scaled by 2^-m.
 * That is, k_mx_exp return z with correction e, such that
 *		exp(x) = (z+e)*2^m
 *
 * Let nb be the number of bit of a table (i.e., sizeof(table)=2^nb).
 *	1. Argument Reduction: given the input x, find r and integer k
 *	   and j such that
 *             	x = (k+j/ts)*(ln2) + r,  |r| <= (1/(2*ts))*ln2 .
 *	   Here ts = 1<<nb is the size of table.
 *
 *	   Method:
 *		n= x*(ts/ln2)+copysign(0.5,x), k = n>>nb; j=n&((1<<nb)-1).
 *		t = x-n*(ln2/ts) .. need multi-precision arithmetic
 *
 * How many extra bit nb is needed in order to guarantee 77 bits accuracy
 * of (24 more bits) r? Since n<=2^(10+nb), we try all n*ln2 (with nb<=7) and
 * found that when n = 116877, n*ln2 is the closest one to a 53 bit floating
 * point number (only the significant parts are displayed)
 *	     116877*ln2  1.B3C74F688A1382000033BCA50B1C3
 *	     rounded to	 1.B3C74F688A1382
 *					 ^^^^^ 18 binary zeros
 *
 * Therefore, to guarantee 77 bits r, we need ln2/ts up to 77+53+18+1 = 149 bits
 * Since n is 10+nb bit <= 17 bit, we break ln2/ts into chunk of 36 bits to
 * guarantee the exactness of the multiplication with n (except the last one).
 * Thus four chunks is needed :36+36+36+53=161 (note that the last one is a
 * full 53 bits precision fp)
 *		t1 = x - n*ln2_nbh1	...exact
 *		t2 = n*ln2_nbh2
 *		t3 = n*ln2_nbh3
 *		t4 = n*ln2_nbh4
 *		ph = t1-t2
 *		pt = t2-(t1-ph)		... ph-pt = t1-t2
 *		if (pt==0.0)	... t1-t2 is exact
 *			qh = ph - t3
 *			qt = t3 - (ph - qh)
 *			if(qt==0.0)	... ph-t3 is exact
 *				r  = qh-t4
 *				rh = (double)(float)r
 *				rt = -(t4-(qh-rh))
 *			else
 *				qt = -(qt+t4);
 *				r = qh+qt;
 *				rh = (double)(float)(r)
 *				rt = qt-(rh-qh)
 *			endif
 *		else
 *			pt  = t3-(pt-t4)
 *			r   = ph+pt
 *			rh  = (double)(float)r
 *			rt  = pt-(rh-ph)
 *		endif
 *
 *	2. exp(x) = 2^k * (2^(j/ts)*exp(r))
 *	   Note:
 *	   a. exp(r) is computed by talyor serier
 *		1+r*(1+r*(1/2!+r*(1/3!+r*(1/4!+...r*(1/10!+....)))))
 *					^double precision arithmetic use here
 *
 * Given nb, how many terms are needed to guarantee 77 bits accuracy.  And also
 * from what term on can double precision arithmetic be used? It depends on the
 * size of (r^k/k!). If r^(k+1)/(k+1)! < 2^-77, then can ignore the k+1 term
 * in summing the taylor serie. And also if r^k/k! <= 2^(-24), then we can
 * apply double precision arithmetic for summing the terms after the k-th terms.
 *		nb	min|r^k/k!|<= 2^-24	min|r^(k+1)/(k+1)!|<2^-77
 *		1	6			14	(remainder: -78.2)
 *		2	5			12	(remainder: -78.4)
 *		3	5			11	(remainder: -83.2)
 *	-->	4	4			9	(remainder: -77.1)
 *		5	4			8	(remainder: -77.2)
 *	-->	6	3			8	(remainder: -86.2)
 *		7	3			7	(remainder: -83.5)
 *		8	3			6	(remainder: -79.0)
 *
 * So from above, only two situation worth to consider: nb=4 and nb=6.
 *
 *	   b. 2^(j/ts) is represented as S1[j]+S2[j], where
 *	      S1[j] = 2^(j/ts) chopped to 24 bits and
 *	      S2[j] = 2^(j/ts) - S1[j].
 *
 * Special cases:
 *	exp(INF) is INF, exp(NaN) is NaN;
 *	exp(-INF)=  0;
 *	for finite argument, only exp(0)=1 is exact.
 *
 * Accuracy:
 *	according to an error analysis, the error is always somewhat less than
 *	2^-22 ulp (unit in the last place).
 *
 * Misc. info.
 *	For IEEE double
 *		if x >  7.09782712893383973096e+02 then exp(x) overflow
 *		if x < -7.45133219101941108420e+02 then exp(x) underflow
 *
 * Constants:
 * The hexadecimal values are the intended ones for the following constants.
 * The decimal values may be used, provided that the compiler will convert
 * from decimal to binary accurately enough to produce the hexadecimal values
 * shown.
 */

#include "mcr.h"

#define	NB 4

#ifdef DEBUG
#define	HA(x)  *(int *)&x, *(1+(int *)&x), *(2+(int *)&x), *(3+(int *)&x)
#include <stdio.h>
#endif

static const double
half[] = {0.5, -0.5},
/* default nb=4 */
invln2_nb = 2.30831206542234141921e+01,	/* 40371547 652B82FE */
ln2_nbh1 = 4.33216951787471771240e-02,	/* 3FA62E42 E0000000 */
ln2_nbh2 = 3.60624929918174075283e-09,	/* 3E2EFA39 E0000000 */
ln2_nbh3 = 1.05532807662813641605e-16,	/* 3C9E6AF2 60000000 */
ln2_nbh4 = 5.15444872131933194148e-24,	/* 3B18ECE6 00FCBDAC */
/* Redundant exception handling code: see __libmcr_exp */
#if 0
o_threshold = 7.09782712893383973096e+02,	/* 0x40862E42, 0xFEFA39EF */
u_threshold = -7.45133219101941108420e+02,	/* 0x40874910, 0xD52D3051 */
twom1000 = 9.33263618503218878990e-302,	/* 0x01700000, 0x00000000 */
huge = 1.0e300,
#endif
one = 1.0,
zero = 0.0;

/*
 * tailor series approximation on [-(ln2)/32,(ln2)/32] exp(x) = 1 + x +
 * 0.5*x*x + ... = 1 + x*( 1 + x*(1/2 + x*(1/6 + x*(1/24 + ...)))) =
 * F0+x*(F1+x*(F2+x*(F3+...)))
 */
#define	F4 ff[8]
#define	F5 ff[9]
#define	F6 ff[10]
#define	F7 ff[11]
#define	F8 ff[12]
#define	F9 ff[13]
static double ff[] = {
	1.0,
	0.0,
	1.0,
	0.0,
	0.5,
	0.0,
	1.66666666666666657415e-01,	/* f[6]=1/6 head 3FC55555 55555555 */
	9.25185853854297065662e-18,	/* f[7]=1/6 tail 3C655555 55555555 */
	4.16666666666666643537e-02,	/* f[8]=1/24  3FA55555 55555555 */
	8.33333333333333321769e-03,	/* 3F811111 11111111 */
	1.38888888888888894189e-03,	/* 3F56C16C 16C16C17 */
	1.98412698412698412526e-04,	/* 3F2A01A0 1A01A01A */
	2.48015873015873015658e-05,	/* 3EFA01A0 1A01A01A */
	2.75573192239858925110e-06,	/* 3EC71DE3 A556C734 */
};

/* nb = 4 */
static const double S1[] = {
	1.00000000000000000000e+00,	/* 3FF00000 00000000 */
	1.04427373409271240234e+00,	/* 3FF0B558 60000000 */
	1.09050762653350830078e+00,	/* 3FF172B8 20000000 */
	1.13878858089447021484e+00,	/* 3FF2387A 60000000 */
	1.18920707702636718750e+00,	/* 3FF306FE 00000000 */
	1.24185776710510253906e+00,	/* 3FF3DEA6 40000000 */
	1.29683947563171386719e+00,	/* 3FF4BFDA C0000000 */
	1.35425543785095214844e+00,	/* 3FF5AB07 C0000000 */
	1.41421353816986083984e+00,	/* 3FF6A09E 60000000 */
	1.47682607173919677734e+00,	/* 3FF7A114 60000000 */
	1.54221081733703613281e+00,	/* 3FF8ACE5 40000000 */
	1.61049032211303710938e+00,	/* 3FF9C491 80000000 */
	1.68179273605346679688e+00,	/* 3FFAE89F 80000000 */
	1.75625205039978027344e+00,	/* 3FFC199B C0000000 */
	1.83400797843933105469e+00,	/* 3FFD5818 C0000000 */
	1.91520655155181884766e+00,	/* 3FFEA4AF A0000000 */
};
static const double S2[] = {
	0.00000000000000000000e+00,	/* 00000000 00000000 */
	4.83347014379782142328e-08,	/* 3E69F312 1EC53172 */
	1.06131749358425753402e-07,	/* 3E7C7D51 7ADCDF7C */
	5.38622214388600755736e-08,	/* 3E6CEAC4 70CD83F5 */
	3.79763538792174980894e-08,	/* 3E64636E 2A5BD1AB */
	4.49683815095311756138e-08,	/* 3E682468 446B6824 */
	7.90192957987462450228e-08,	/* 3E75362A 271D4397 */
	1.09085940579860513649e-07,	/* 3E7D4854 2958C930 */
	2.42032342089579361800e-08,	/* 3E59FCEF 32422CBE */
	7.42003025340431556211e-08,	/* 3E73EB01 86D7D510 */
	8.07090469079979051284e-09,	/* 3E415506 DADD3E2A */
	9.83621719880451981717e-09,	/* 3E451F84 80E3E235 */
	9.44539622891872401233e-08,	/* 3E795AD3 AD5E8734 */
	1.09973519209674659656e-07,	/* 3E7D8552 9C2220CB */
	1.07970011408799574157e-07,	/* 3E7CFBA4 8725DA05 */
	9.84532844621636118964e-09,	/* 3E452486 CC2C7B9D */
};

double
__libmcr_k_mx_exp(x, e, m)
	double x, *e;
	int *m;			/* exp(x) = 2^m * m2_exp() */
{
	double y, c, t, t1, t2, t3, t4, p;
	double qh, qt, ph, pt, rh, rt, r, ee[2], z[3];

#ifdef DEBUG
	long double lx, ly, lz, expl(), logl(), lw, lp, exp2l();
	double w;
	int i;
#endif
	int k, xsb, n, j;
	int hx;

	hx = HIGH_WORD(x);	/* high word of x */
	xsb = (hx >> 31) & 1;	/* sign bit of x */
	hx = hx & 0x7fffffff;	/* high word of |x| */
	*e = 0.0;
	*m = 0;

/* Redundant exception handling code: see __libmcr_exp */
#if 0
	/* filter out non-finite argument */
	if (hx >= 0x40862E42) {	/* if |x|>=709.78... */
		if (hx >= 0x7ff00000) {
			if (((hx & 0xfffff) | LOW_WORD(x)) != 0)
				return (x + x);	/* NaN */
			else
			{
				/* exp(+-inf)={inf,0} */
				return ((xsb == 0) ? x : 0.0);
			}
		}
		if (x > o_threshold)
			return (huge * huge);	/* overflow */
		if (x < u_threshold)
			return (twom1000 * twom1000);	/* underflow */
	}
#endif

	if (hx <= 0x3e600000) {	/* |x|<= 2^-25 */
		if (hx <= 0x3d900000) {	/* |x|<=2^-38 */
			t = one + x;
			*e = x - (t - one);
		} else {
			t = one + x * (one + x * 0.5);
			*e = x * x * 0.5 - ((t - one) - x);
		}
		return (t);
	}
	/* compute (rh+rt) = x - n*(ln2)/16 */
	n = x * invln2_nb + half[xsb];	/* determine n by rounding x*16/ln2 */
	k = n >> NB;
	j = n & ((1 << NB) - 1);	/* k = [n/16], j = n-16k */
	c = (double) n;
	/* simulate higher precision arithmetic */
	t1 = x - c * ln2_nbh1;
	t2 = c * ln2_nbh2;
	ph = t1 - t2;
	t3 = c * ln2_nbh3;
	pt = t2 - (t1 - ph);
	t4 = c * ln2_nbh4;
	if (pt == zero) {
		qh = ph - t3;
		qt = t3 - (ph - qh);
		ph = qh;
		if (qt == zero)
			pt = -t4;
		else
			pt = -(qt + t4);
	} else
		pt = t3 - (pt - t4);
	/* represent result in two floating point numbers */
	r = ph + pt;
	rh = (double) (float) r;
	rt = pt - (rh - ph);
	/* fast computation on the tail part of the polynomial */
	ee[0] = r * (F4 + r * (F5 + r * (F6 + r * (F7 + r * (F8 + r * F9)))));
	z[0] = r;
	z[1] = rh;
	z[2] = rt;
#ifdef DEBUG
	i = 1 << NB;
	w = 1.0 / (double) i;
	lx = (long double) x - logl(2.0L) * (long double) (w * (double) n);
	ly = (long double) rh + (long double) rt;
	printf(" r    = %08X%08X\n", r);
	printf("rh+rt =%08X%08X %08X %08X\n", HA(ly));
	printf("lr    =%08X%08X %08X %08X\n", HA(lx));
	lz = expl(lx);
	lw = (lz - 1.0L - lx - lx * lx * 0.5L) / (lx * lx);
	lw = (lw - lx / 6.0L) / lx;
#endif
	__libmcr_mx_poly(z, ff, ee, 3);

	/* 2^(j/ts)*(ph+pt) */
	t = S1[j] * ee[0];
	p = S2[j] * ee[0] + (S1[j] + S2[j]) * ee[1];
	y = p + t;
	*e = p - (y - t);
#ifdef DEBUG
	printf("   S1[%2d] = %08X %08X\n", j, S1[j]);
	printf("   S2[%2d] = %08X %08X\n", j, S2[j]);
	ly = (long double) S1[j] + (long double) S2[j];
	printf("     S1+S2 =%08X%08X %08X %08X\n", HA(ly));
	lp = exp2l((long double) ((double) j / 64.0));
	printf("exp2l(j/64)=%08X%08X %08X %08X\n", HA(lp));
	lw *= lp;
	ly = (long double) y + (long double) *e;
	printf("final exp (unscale)  =%08X%08X %08X %08X\n", HA(ly));
	printf("longd exp (unscale)  =%08X%08X %08X %08X\n", HA(lw));
#endif
	*m = k;
	return (y);
}
