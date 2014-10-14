
#pragma ident "@(#)__libmcr_mx_atan.c 1.6 04/02/25 SMI"

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
 * __libmcr_mx_atan(x,err)
 * Table look-up algorithm
 * By K.C. Ng, March 9, 1989
 *
 * Algorithm.
 *
 * The algorithm is based on atan(x)=atan(y)+atan((x-y)/(1+x*y)).
 * We use poly1(x) to approximate atan(x) for x in [0,1/8] with
 * error (relative)
 * 	|(atan(x)-poly1(x))/x|<= 2^-83.41
 *
 * and use poly2(x) to approximate atan(x) for x in [0,1/65] with
 * error
 *	|atan(x)-poly2(x)|<= 2^-86.8
 *
 * Here poly1 and poly2 are odd polynomial with the following form:
 *		x + x^3*(a1+x^2*(a2+...))
 *
 * (0). Purge off Inf and NaN and 0
 * (1). Reduce x to positive by atan(x) = -atan(-x).
 * (2). For x <= 1/8, use
 *	(2.1) if x < 2^(-prec/2), atan(x) =  x  with inexact flag raised
 *	(2.2) Otherwise
 *		atan(x) = poly1(x)
 * (3). For x >= 8 then (prec = 78)
 *	(3.1) if x >= 2^prec,     atan(x) = atan(inf) - pio2lo
 *	(3.2) if x >= 2^(prec/3), atan(x) = atan(inf) - 1/x
 *	(3.3) if x >  65,         atan(x) = atan(inf) - poly2(1/x)
 *	(3.4) Otherwise,	  atan(x) = atan(inf) - poly1(1/x)
 *
 * (4). Now x is in (0.125, 8)
 *      Find y that match x to 4.5 bit after binary (easy).
 *	If iy is the high word of y, then
 *		single : j = (iy - 0x3e000000) >> 19
 *		double : j = (iy - 0x3fc00000) >> 16
 *		quad   : j = (iy - 0x3ffc0000) >> 12
 *
 *	Let s = (x-y)/(1+x*y). Then
 *	atan(x) = atan(y) + poly1(s)
 *		= _TBL_atan_hi[j] + (_TBL_atan_lo[j] + poly2(s) )
 *
 *	Note. |s| <= 1.5384615385e-02 = 1/65. Maxium occurs at x = 1.03125
 *
 */

#include "mcr.h"

#define	P1 p[2]
#define	P4 p[8]
#define	P5 p[9]
#define	P6 p[10]
#define	P7 p[11]
#define	P8 p[12]
#define	P9 p[13]
static double p[] = {
	1.0,
	0.0,
	-3.33333333333333314830e-01,	/* p1   = BFD55555 55555555 */
	-1.85030852238476921863e-17,	/* p1_l = BC755525 9783A49C */
	2.00000000000000011102e-01,	/* p2   = 3FC99999 9999999A */
	-1.27263196576150347368e-17,	/* p2_l = BC6D584B 0D874007 */
	-1.42857142857141405923e-01,	/* p3   = BFC24924 9249245E */
	-1.34258204847170493327e-17,	/* p3_l = BC6EF534 A112500D */
	1.11111111110486909803e-01,	/* p4   = 3FBC71C7 1C71176A */
	-9.09090907557387889470e-02,	/* p5   = BFB745D1 73B47A7D */
	7.69230541541713053189e-02,	/* p6   = 3FB3B13A B1E68DE6 */
	-6.66645815401964159097e-02,	/* p7   = BFB110EE 1584446A */
	5.87081768778560317279e-02,	/* p8   = 3FAE0EFF 87657733 */
	-4.90818147456113240690e-02,	/* p9   = BFA92140 6A524B5C */
};
#define	Q1 q[2]
#define	Q3 q[6]
#define	Q4 q[7]
#define	Q5 q[8]
static double q[] = {
	1.0,
	0.0,
	-3.33333333333333314830e-01,	/* q1   = BFD55555 55555555 */
	-1.85022941571278638733e-17,	/* q1_l = BC7554E9 D20EFA66 */
	1.99999999999999927836e-01,	/* q2   = 3FC99999 99999997 */
	-1.28782564407438833398e-17,	/* q2_l = BC6DB1FB 17217417 */
	-1.42857142855492280642e-01,	/* q3   = BFC24924 92483C46 */
	1.11111097130183356096e-01,	/* q4   = 3FBC71C6 E06595CC */
	-9.08553303569109294013e-02,	/* q5   = BFB7424B 808CDA76 */
};
static const double
one = 1.0,
pio2hi = 1.570796326794896558e+00,
pio2lo = 6.123233995736765886e-17;

double
__libmcr_mx_atan(x, err)
	double x, *err;
{
	double y, z, r, s, t, w, s_h, s_l, x_h, x_l, zz[3], ee[2], z_h,
		z_l, r_h, r_l, u, v;
	int ix, iy, sign, j;

	ix = HIGH_WORD(x);
	sign = ix & 0x80000000;
	ix ^= sign;

	/* for |x| < 1/8 */
	if (ix < 0x3fc00000) {
		if (ix < 0x3f300000) {	/* when |x| < 2**-12 */
			if (ix < 0x3d800000) {	/* if |x| < 2**-39 */
				*err = (double) ((int) x);
				return (x);
			}
			z = x * x;
			t = x * z * (q[2] + z * (q[4] + z * q[6]));
			r = x + t;
			*err = t - (r - x);
			return (r);
		}
		z = x * x;

		/* use double precision at p4 and on */
		ee[0] = z *
			(P4 + z *
			(P5 + z * (P6 + z * (P7 + z * (P8 + z * P9)))));

		x_h = (double) ((float) x);
		z_h = (double) ((float) z);
		x_l = x - x_h;
		z_l = (x_h * x_h - z_h);
		zz[0] = z;
		zz[1] = z_h;
		zz[2] = z_l + x_l * (x + x_h);

		/*
		 * compute (1+z*(p1+z*(p2+z*(p3+e)))) by call
		 * __libmcr_mx_poly
		 */

		__libmcr_mx_poly(zz, p, ee, 3);

		/* finally x*(1+z*(p1+...)) */
		r = x_h * ee[0];
		t = x * ee[1] + x_l * ee[0];
		s = t + r;
		*err = t - (s - r);
		return (s);
	}
	/* for |x| >= 8.0 */
	if (ix >= 0x40200000) {	/* x >=  8 */
		x = fabs(x);
		if (ix >= 0x42600000) {	/* x >=  2**39 */
			if (ix >= 0x44c00000) {	/* x >=  2**77 */
/* Redundant exception handling code: see __libmcr_atan */
#if 0
				if (ix >= 0x7ff00000) {
				    if (((ix-0x7ff00000) | LOW_WORD(x)) != 0) {
					*err = 0.0;
					return (x - x);
				    }
				}
#endif
				y = -pio2lo;
			} else
				y = one / x - pio2lo;
			if (sign == 0) {
				t = pio2hi - y;
				*err = -(y - (pio2hi - t));
			} else {
				t = y - pio2hi;
				*err = y - (pio2hi + t);
			}
			return (t);
		} else {
			/* compute r = 1/x */
			r = one / x;
			z = r * r;
			if (ix < 0x40504000) {	/* 8 <  x <  65 */

				/* use double precision at p4 and on */
				ee[0] = z *
					(P4 + z *
					(P5 + z *
					(P6 + z * (P7 + z * (P8 + z * P9)))));
				x_h = (double) ((float) x);
				r_h = (double) ((float) r);
				z_h = (double) ((float) z);
				r_l = r * ((x_h - x) * r_h - (x_h * r_h - one));
				z_l = (r_h * r_h - z_h);
				zz[0] = z;
				zz[1] = z_h;
				zz[2] = z_l + r_l * (r + r_h);
				/*
				 * compute (1+z*(p1+z*(p2+z*(p3+e)))) by call
				 * __libmcr_mx_poly
				 */
				__libmcr_mx_poly(zz, p, ee, 3);
			} else { /* x < 65 < 2**39 */
				/* use double precision at q3 and on */
				ee[0] = z * (Q3 + z * (Q4 + z * Q5));
				x_h = (double) ((float) x);
				r_h = (double) ((float) r);
				z_h = (double) ((float) z);
				r_l = r * ((x_h - x) * r_h - (x_h * r_h - one));
				z_l = (r_h * r_h - z_h);
				zz[0] = z;
				zz[1] = z_h;
				zz[2] = z_l + r_l * (r + r_h);
				/*
				 * compute (1+z*(q1+z*(q2+e))) by call
				 * __libmcr_mx_poly
				 */
				__libmcr_mx_poly(zz, q, ee, 2);
			}
			/* pio2 - r*(1+...) */
			v = r_h * ee[0];
			t = pio2lo - (r * ee[1] + r_l * ee[0]);
			if (sign == 0) {
				s = pio2hi - v;
				t -= (v - (pio2hi - s));
			} else {
				s = v - pio2hi;
				t = -(t - (v - (s + pio2hi)));
			}
			w = s + t;
			*err = t - (w - s);
			return (w);
		}
	}
	/* now x is between 1/8 and 8 */
	HIGH_WORD(x) = ix;
	iy = (ix + 0x00008000) & 0x7fff0000;
	HIGH_WORD(y) = iy;
	LOW_WORD(y) = 0;
	j = (iy - 0x3fc00000) >> 16;

	w = (x - y);
	v = 1 / (one + x * y);
	s = w * v;
	z = s * s;
	/* use double precision at q3 and on */
	ee[0] = z * (Q3 + z * (Q4 + z * Q5));
	s_h = (double) ((float) s);
	z_h = (double) ((float) z);
	x_h = (double) ((float) x);
	t = (double) ((float) (one + x * y));
	r = -((x_h - x) * y - (x_h * y - (t - one)));
	s_l = -v * (s_h * r - (w - s_h * t));
	z_l = (s_h * s_h - z_h);
	zz[0] = z;
	zz[1] = z_h;
	zz[2] = z_l + s_l * (s + s_h);
	/* compute (1+z*(q1+z*(q2+e))) by call __libmcr_mx_poly */
	__libmcr_mx_poly(zz, q, ee, 2);
	v = s_h * ee[0];
	t = __libmcr_TBL_atan_lo[j] + (s * ee[1] + s_l * ee[0]);
	u = __libmcr_TBL_atan_hi[j];
	s = u + v;
	t += (v - (s - u));
	w = s + t;
	*err = t - (w - s);
	if (sign != 0) {
		w = -w;
		*err = -*err;
	}
	return (w);
}
