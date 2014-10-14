
#pragma ident "@(#)__libmcr_mx_log.c 1.5 04/02/25 SMI"

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
 * double __libmcr_mx_log((double)x,(double*)e)
 * Modified multiprecision log with error correction e up to 22 extra bits.
 * Method:
 *	(Assume x is normalize and  > 0)
 *	1. Argument Reduction. Find i,n,and r such that
 *		x = 2**n r,
 *	   where HI(r) 	is between 0x3fe78000 and 0x3ff77fff
 *		 i 	is the index such that (0x3fe8+i)<<16 is close
 *			to HI(r).
 *	   Implementation note:
 *		j    = (hx>>15)+0x11;
 *		n    = (j>>5)-0x3ff;
 *		hx  -= (n<<20);
 *		i    = ((hx+0x8000)>>16)-0x3fe8;
 *
 *	   Let y[i] be a number whose leading word is close to
 *	   (0x3fe8+i)<<16 and be such that log(y[i]) is almost
 *	   in double precision (with 24 extra 0 bits). Also,
 *	   Max(|y-r|)  <= 2**-5 + 2**-12 < 0.0315.
 *
 *	   We have
 *		log(x) = n*log2 + log(y[i]) + log(r/y[i])
 *
 *	2. Taylor series for log(r/y[i]).
 *	   Let y stand for y[i] from now on.
 *	   Let s = (r-y)/(r+y). Then
 *		log(r/y) = log((1+s)/(1-s)) = 2s + 2/3 s^3 + 2/5 s^5 + ...
 *	   Since |r-y|<=0.0315, we have |s|<0.01575, s*s <0.0002480625
 *	   = 2**(-11.977..).
 *	   7 terms is enough for 77 bits accuracy, and one can use double
 *	   precision arithemtic from the third term on (the third term has
 *	   more than 24 bits shift to the left of s.)
 *
 *			  3        5        7        9         11        13
 *		2s + 2/3 s  + 2/5 s  + 2/7 s  + 2/9 s  + 2/11 s  + 2/13 s
 *			       ^^^^^^^^^^^^^^^^
 *			      from this terms on use double precision arithmetic
 *	   Note:
 *		How to compute s1 + s2 = s, where s1 is the first 24 bits of s.
 *		Let p = r-y (exact), q = r+y (rounded), q1 = 27 bits of q,
 *		and q2 = (r+y)-q1 = (r-q1)+y; (q2 has at most 27 bits)
 *		Let qi = 1/q. Then
 *			s = p*qi; s1 = s rounded to 24 bit;
 *			s2 = p/(q1+q2)-s1 = qi*(p-s1*q1-s1*q2)
 *			   = qi*((p-s1*q1)-s1*q2)
 *		How to compute z1+z2= (s1+s2)**2
 *			z1 = s*s rounded to 24 bits
 *			z2 = (s1*s1+2s1*s2+s2*s2)-z1
 *			   = ((s1*s1-z1)+2s1*s2)+s2*s2:
 *
 *		Then log(r/y) is computed by
 *			(s1+s2)*(2.0+(z1+z2)*(2/3+....)) = g1 + g2
 *
 *	3. Table look up for log(y[i])
 *		j    = (hx>>15)+0x11;
 *		n    = (j>>5)-0x3ff;
 *		hx  -= (n<<20);
 *		i    = ((hx+0x8000)>>16)-0x3fe8;
 *		y      = ylny[2i]
 *		log(y) = ylny[2i+1]
 *
 *	4. n*ln(2)
 *	   since n at most 11 bits, let ln2_h = ln2 chopped to 42 bits,
 *	   ln2_l = (ln2-ln2_h). Then we have ~42+53 = 95 accuracy.
 *
 *
 *	5. Sum them together
 *		[      n*ln_h      ] [      n*ln_l      ]
 *		     [     ln(y  )    ]
 *			 [   g1          ] [     g2            ]
 *		--------------------------------------------------
 *	   (i) 	A = n*ln_h + ln(y) + g1
 *		B = n*ln_l + g2 + (((n*ln_h-A)+ln(y))+g1)
 *		result = A+B;
 *		error = B-(result-A)
 *
 */

#include "mcr.h"

static const double ylny[] = {	/* y = ylny[2i], log(y) = ylny[2i+1] */
	7.50085672710530171337e-01, -2.87567848694855743297e-01,
	7.81156543355310373222e-01, -2.46979709592310303634e-01,
	8.12438483237547193205e-01, -2.07715080583017158711e-01,
	8.43794452953426921127e-01, -1.69846353201363070573e-01,
	8.75109166773791691263e-01, -1.33406638379466563338e-01,
	9.06227439459201256078e-01, -9.84649675129698753739e-02,
	9.37449804002442932394e-01, -6.45920649684107034405e-02,
	9.68854432006046706327e-01, -3.16409033442315718032e-02,
	1.0e+00, 0.0,
	1.06248367360119333469e+00, 6.06092556759700376579e-02,
	1.12498771534588781762e+00, 1.17772115904219015770e-01,
	1.18751408372130939917e+00, 1.71862116832169758984e-01,
	1.24985120463076171404e+00, 2.23024507933437032836e-01,
	1.31270013179742672804e+00, 2.72086185229007349040e-01,
	1.37507504991287965446e+00, 3.18508311383826403507e-01,
	1.43752755519066410805e+00, 3.62924662333936887126e-01,
};
static const double
ln2_h = 6.93147180559890330187e-01,	/* 3FE62E42 FEFA3800 */
ln2_l = 5.49792301870837115524e-14,	/* 3D2EF357 93C76730 */
two54 = 18014398509481984.0,
/* Redundant exception handling code: see __libmcr_log */
#if 0
one = 1.0,
zero = 0.0,
#endif
v1 = 2.0 / 5.0,
v2 = 2.0 / 7.0,
v3 = 2.0 / 9.0,
v4 = 2.0 / 11.0,
v5 = 2.0 / 13.0,
c2_3_h = 2.0 / 3.0,
c2_3_l = ((2.0 - 4.0 / 3.0) - 2.0 / 3.0) / 3.0;

double
__libmcr_mx_log(double x, double *err)
{
	int hx, i, j, n, na;
	double y, p, q, q1, qi, q2, s, f1, f2, u, w, s1, s2, w1, w2, z1,
		z2, t1, t2, g1, g2, t;

	hx = HIGH_WORD(x);
	*err = 0.0;
	na = 0;
	/* filter out exception argument */
/* Redundant exception handling code: see __libmcr_log */
#if 0
	if ((hx + 0x00208000) < 0x00308000) {
		/* subnormal,0,negative,inf,nan,huge */
		i = ((unsigned) (hx << 1)) >> 1;
		if (hx >= 0x7ff00000)
			return (x + x);	/* x is inf/nan */
		if ((i | LOW_WORD(x)) == 0)
			return (-one / (x * x));	/* log(+-0) is -inf */
		if (hx < 0)
			return ((x - x) / zero);	/* log(-#) is NaN */
#endif
		if (hx < 0x00100000) {	/* subnormal */
			x *= two54;
			hx = HIGH_WORD(x);
			na = -54;
		}
/* Redundant exception handling code: see __libmcr_log */
#if 0
	}
#endif
	j = (hx >> 15) + 0x11;
	n = (j >> 5) - 0x3ff;
	hx -= (n << 20);
	HIGH_WORD(x) = hx;
	i = ((hx + 0x8000) >> 16) - 0x3fe8;
	na += n;
	f1 = na * ln2_h;
	f2 = na * ln2_l;
	i = i + i;
	y = ylny[i];
	/* log(x/y) */
	q = x + y;
	p = x - y;
	q1 = q;
	LOW_WORD(q1) &= 0xfc000000;
	qi = 1.0 / q;
	q2 = y - (q1 - x);
	s = p * qi;
	u = s * s;
	w = c2_3_l + u * (v1 + u * (v2 + u * (v3 + u * (v4 + u * v5))));
	s1 = (double) ((float) s);
	s2 = -qi * (s1 * q2 - (p - s1 * q1));
	w1 = (double) ((float) (w + c2_3_h));
	w2 = w - (w1 - c2_3_h);
	z1 = (double) ((float) u);
	z2 = (((s1 + s1) * s2) - (z1 - s1 * s1)) + s2 * s2;
	/*
	 * log(x/y) = g1+g2 = (s1+s2)*(2.0+(z1+z2)*(w1+w2))
	 */
	t = z1 * w1;
	w = w1 * z2 + u * w2;
	t1 = (double) ((float) (t + w + 2.0));
	t2 = (t - (t1 - 2.0)) + w;
	p = s1 * t1;
	q = s2 * t1 + s * t2;
	g1 = p + q;
	g2 = q - (g1 - p);
	t1 = ylny[i + 1];
	z1 = f1 + t1 + g1;
	f2 += g2 - (((z1 - f1) - t1) - g1);
	s = z1 + f2;
	*err = f2 - (s - z1);
	return (s);
}
