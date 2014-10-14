
#pragma ident "@(#)__libmcr_pow.c 1.10 04/04/09 SMI"

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

#ifdef DEBUG
#include <stdio.h>
#endif
/*
 * correctly rounded pow
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
 * Additional special cases for double precision arithmetic
 *	1.   x ** -1    return  1/x
 *	2.   x ** 2     return  x*x
 *	3.  +x ** 1/2   return  sqrt(x)
 *	4.  if y > 2**64, x**y  overflow if x>1, underflow if x<1
 *	5.  if y > 2**31, x**y  overflow  if Hi(x)>=0x3ff00001;
 *				underflow if Hi(x)<=0x3feffffe;
 *	6.  Exact and half-way cases of x**y
 *		(No more than 54 bits mantissa result)
 *
 *	    Let x = Nx * 2**(kx), y = Ny * 2**(ky), and z = Nz * 2**(kz), where
 *	    Nx,Ny,Nz are odd integer and Nx,Ny < 2**53, Nz < 2**54.
 *
 *	    If z = x**y, we have
 *		Nz * 2**(kz) = [Nx * 2**(kx)]**(Ny*2**ky)
 *			     = Nx**(Ny*2**ky) * 2**(kx*Ny*2**ky)
 *		==>
 *			(A)	Nz = Nx**(Ny*2**ky) and
 *			(B)	kz = kx*Ny*2**ky.
 *	    (1) Nx = 1.
 *		We have x = 2**kx and z = 2**(kx*y). One needs only to check
 *		whether kx*y is an integer within bound.
 *	    (2) Nx = p**px for some prime number p>2 and positive integer px.
 *		Let us consider the cases when ky < 0.
 *		In this case, since kz is an integer and Ny is odd, then (B)
 *		implies 2**(-ky) | kx and (A) implies Nx**(2**(ky)) is an odd
 *		integer j.  Thus x = Nx * 2**kx can be square rooted -ky
 *		times exactly. How large could ky be? Since Nx = j**(2**(-ky))
 *		for some odd integer j>1, we have Nx >= 3**(2**(-ky)).
 *		Since Nx < 2**53 < 3**(2**(6)), (-ky) <= 5.
 *		(if x is subnormal, (-ky) could go up to 1024.)
 *
 *		Thus, we can apply the following algorithm to reduce negative
 *		ky to 0:
 *		    If y is not an integer and 1024*y is an integer, then do:
 *			1. Determine k>0 such that 2**k * y is an odd integer.
 *			2. Clear inexact bit and let w = x.
 *			3. Perform w = sqrt(w) k time
 *			4. check if all sqrt operations are exact. If yes,
 *			   replace x by w and y by y/2**k. Otherwise, don't
 *			   do anything.
 *			5. Continue.
 *		One can also check every sqrt() and stop when inexact bit is on.
 *
 *		If the access of inexact bit is not available, then one must
 *		check whether the squaring of w is exactly equal to the previous
 *		one.
 *
 *		Let z be the previous number and w = sqrt(z),
 *		and let nz = ilogb(z)
 *		and nw = ilogb(w). If mw deonote the length of non-zero mantissa
 *		of z, then
 *			w is the exact sqrt of z
 *		if and only if
 *			(mw <=26 or (mw=27 and nw+nw=nz)) and w*w=z
 *
 *	    (3) Now we can concentrate on y is an integer and Nx is >= 3. When
 *		y is 2, pow just return x*x. Since Nx**34 >= 3**34 has exactly
 *		54 significant bits, we need only to concern 3 <= y <= 34.
 *
 *		For each particular integer y that is between 3 and 34, the size
 *		of Nx is limited by (2**54)**(1/y).
 *		Below is a tabulated table for the upper bound of
 *		Nx for various y:
 *
 *			y	upper bound of Nx
 *		     -----------------------------
 *			3	262143
 *			4	11585
 *			5	1781
 *			6	511
 *			7	209
 *			8	107
 *			9	63
 *			10	41
 *			11	29
 *			12	21
 *			13	17
 *			14	13
 *			15	11
 *			16,17	9
 *			18-19	7
 *			20-23	5
 *			24-34	3
 *
 *	    (4) Compute x**y by repeat multiplication. (That way will preserve
 *		exactness and rounded the half-way case correctly.)
 *
 * Constants :
 * The hexadecimal values are the intended ones for the following
 * constants. The decimal values may be used, provided that the
 * compiler will convert from decimal to binary accurately enough
 * to produce the hexadecimal values shown.
 */


#include "mcr.h"

#pragma weak pow = __libmcr_pow

#ifdef DEBUG
#define	H2(x)  *(int *)&x, *(1+(int *)&x)
#endif

#ifdef DEBUG
extern int ieee_flags();
#endif

static const double
rndc = -1.02735137937997933477e+00,	/* 0xbff07007 ffff0000 */
huge = 1.0e300,
two1000 = 1.07150860718626732095e+301,	/* 0x7E700000 00000000 */
two100 = 1.26765060022822940150e+30,	/* 0x46300000 00000000 */
zero = 0.0,
one = 1.0,
tiny = 2.22507385850720187716e-308;	/* 0x00100000 0x00000001 */

double
__libmcr_pow(double x, double y)
{
	int hx, hy, ix, iy, rnd = 0, nw, ly, lx, ncrd, m, i, iz;
	int yisint, k, j, hz;
	int inx, nx;
#ifdef DEBUG
	int in[10], *out;
	int status;
#endif
	double z, er, w, t, ax, s;
	float fw, ft;

	hx = HIGH_WORD(x);
	lx = LOW_WORD(x);
	hy = HIGH_WORD(y);
	ly = LOW_WORD(y);
	ix = hx & 0x7fffffff;
	iy = hy & 0x7fffffff;

#ifdef DEBUG
	printf("Input x,y = %1.20e %1.20e\n", x, y);
	printf("In hex: %08X %08X %08X %08X\n", hx, lx, hy, ly);
#endif
	/* y==zero: x**0 = 1 */
	if ((iy | ly) == 0)
		return (one);

	/* +-NaN return x+y */
	if (ix > 0x7ff00000 || ((ix == 0x7ff00000) && (lx != 0)) ||
	    iy > 0x7ff00000 || ((iy == 0x7ff00000) && (ly != 0)))
		return (x + y);

	/*
	 * determine if y is an odd int yisint = 0	... y is not an
	 * integer yisint = 1	... y is an odd int yisint = 2	... y is an
	 * even int
	 */
	yisint = 0;
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
#ifdef DEBUG
	status = ieee_flags("get", "exception", in, &out);
	printf("yisint = %d  FLAGS = %d\n", yisint, (status & 0x3ff));
#endif
	/* special value of y */
	if ((ly | (iy << 12)) == 0) {
		if (iy == 0x7ff00000) {	/* y is +-inf */
			if (((ix - 0x3ff00000) | lx) == 0)
				return (y - y);	/* inf**+-1 is NaN */
			else if (ix >= 0x3ff00000) {
				/* (|x|>1)**+-inf =  inf,0 */
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
#ifdef DEBUG
	printf("Pass special value of y\n");
#endif

	ax = fabs(x);
	/* special value of x */
	if (lx == 0) {
#ifdef DEBUG
		printf("x = %08X %08X %1.20e\n", H2(x), x);
		printf("j = %08X\n", j);
#endif
		/*
		 * Tried the 2 lines below to replace the following line but
		 * failed at subnormal number j=((ix<<2)>>22); if
		 * ((j+(j&1))==0) {
		 */
		if ((ix << 12) == 0) {
			j = ix >> 20;
			/*
			 * if(ix==0x7ff00000||ix==0||ix==0x3ff00000){
			 */
			if (j == 0x7ff || j == 0 || j == 0x3ff) {
				z = ax;	/* x is +-0,+-inf,+-1 */
				if (hy < 0)
					z = one / z;	/* z = (1/|x|) */
				if (hx < 0) {
					/*
					 * if(((ix-0x3ff00000)|yisint)==0) {
					 */
					if (((j - 0x3ff) | yisint) == 0) {
						/* (-1)**non-int is NaN */
						z = (z - z) / (z - z);
					} else if (yisint == 1) {
						/* (x<0)**odd = -(|x|**odd) */
						z = -z;
					}
				}
				return (z);
			}
		}
		if ((ly | (hy << 14)) == 0) {
			j = (hy >> 18) - 0x1002;
			if (j == 0)
				return ((x * x) * x);	/* y=3 */
			if (j >= 2 && j <= 4) {
				t = x * x;
				if (j == 2)
					return (t * t);	/* y=4 */
				if ((hx & 0xf) == 0) {
					w = x * t;	/* x*x*x exact */
					if (j == 3)
						return (t * w);	/* y=5 */
					else
						return (w * w);	/* y=6 */
				}
			}
#if 0
			/*
			 * the following code is replaced by the above 9
			 * lines
			 */
			if (hy == 0x40080000)
				return ((x * x) * x);	/* y is 3, x*x exact */
			if (hy == 0x40100000) {
				t = x * x;
				return ((t * t));	/* y is 4, x*x exact */
			}
			if (hy == 0x40140000 && (hx & 0xf) == 0) {
				/* y**5 = y**2 * y**3 */
				t = x * x;
				w = x * t;	/* exact */
				return (t * w);
			}
#endif
		}
	}
#ifdef DEBUG
	printf("Pass special value of x\n");
#endif
	/* (x<0)**(non-int) is NaN */
	if ((((hx >> 31) + 1) | yisint) == 0)
		return ((x - x) / (x - x));

	s = one;		/* s (sign of result -ve**odd) = -1 else = 1 */
	if (hx < 0) {		/* (-ve)**(odd int) */
		x = -x;
		if (yisint == 1)
			s = -one;
	}
	/* |y| is huge */
	if (iy > 0x41e00000) {	/* if |y| > 2**31 */
		z = s * huge;
		t = s * tiny;
		if (iy > 0x43f00000) {	/* if |y| > 2**64, must o/uflow */
			if (ix <= 0x3fefffff)
				return ((hy < 0) ? z * huge : t * tiny);
			if (ix >= 0x3ff00000)
				return ((hy > 0) ? z * huge : t * tiny);
		}
		/* over/underflow if x is not close to one */
		if (ix < 0x3fefffff)
			return ((hy < 0) ? z * huge : t * tiny);
		if (ix > 0x3ff00000)
			return ((hy > 0) ? z * huge : t * tiny);
	}
#ifdef DEBUG
	printf("Pass |y| is huge.\n");
#endif

	z = 1024.0 * y;
	i = (int) z;
	if (yisint == 0 && (double) i == z) {
		/*
		 * determine k such that  2 ** k * y is an odd integer,
		 * then change the problem to
		 *	x <- sqrt ** k (x), y <- 2 ** k * y,
		 * provided that sqrt ** k (x) is exact.
		 */
		j = 1024;
		while ((i & 1) == 0) {
			i >>= 1;
			j >>= 1;
		}
		k = j >> 1;
		z = x;
		inx = 0;
#ifdef DEBUG
		printf("k,inx = %d %d, x,y = %1.20e %1.20e\n", k, inx, x, y);
#endif
		while (k > 0 && inx == 0) {
			z = __libmcr_k_exactsqrt(z, &inx);
			k >>= 1;
		}
		if (inx == 0) {
			x = z;
			ix = HIGH_WORD(z);
			lx = LOW_WORD(z);
			y *= (double) j;
			hy = HIGH_WORD(y);
			iy = hy & 0x7fffffff;
			yisint = 1;
		}
#ifdef DEBUG
		printf("inx = %d, x,y = %1.20e %1.20e\n", inx, x, y);
#endif
	}
	/*
	 * Is the result less than 54 bits? (now only happened when y is an
	 * integer)
	 */
	if (yisint > 0) {
#ifdef DEBUG
		status = ieee_flags("get", "exception", in, &out);
		printf("yisint = %d  FLAGS = %d\n", yisint, (status & 0x3ff));
#endif
		if (((ix & 0xfffff) | lx) == 0) { /* x is +-2**n, y integers */
			k = (ix >> 20) - 0x3ff;
			z = y * (double) k;
			if (z >= 1024.0) {
				t = s * two100;
				return (t * two1000);	/* overflow */
			} else if (z <= -1077) {
				t = s / two100;
				return (t / two1000);	/* underflow */
			} else {
				k = (int) z;
				i = k / 2;
				j = k - i;
				w = t = one;
				HIGH_WORD(t) = ((0x3ff) + j) << 20; /* t=2**j */
				HIGH_WORD(w) = ((0x3ff) + i) << 20; /* w=2**i */
				t *= s;
				return (t * w);
			}
		};
		if (iy == 0x3ff00000) {
			if (hy > 0)
				return (s * x);	/* y =  1 */
			else
				return (s / x);	/* y = -1 */
		} else if (hy == 0x40000000)
			return (x * x);	/* y = 2 */
		else if (lx == 0 && hy > 0 && iy <= 0x40410000) {
			/*
			 * only 3<=y<=34 is consider now, since 3**35 have
			 * more than 54 bits
			 */
			i = (int) y;
			nx = (ix & 0x000fffff) | 0x00100000;
			while ((nx & 1) == 0) {
				nx >>= 1;
			}	/* odd x */
			k = 0;
			switch (i) {
			case 3:
				if (nx <= 262143)
					k = 1;
				break;
			case 4:
				if (nx <= 11585)
					k = 1;
				break;
			case 5:
				if (nx <= 1781)
					k = 1;
				break;
			case 6:
				if (nx <= 511)
					k = 1;
				break;
			case 7:
				if (nx <= 209)
					k = 1;
				break;
			case 8:
				if (nx <= 107)
					k = 1;
				break;
			case 9:
				if (nx <= 63)
					k = 1;
				break;
			case 10:
				if (nx <= 41)
					k = 1;
				break;
			case 11:
				if (nx <= 29)
					k = 1;
				break;
			case 12:
				if (nx <= 21)
					k = 1;
				break;
			case 13:
				if (nx <= 17)
					k = 1;
				break;
			case 14:
				if (nx <= 13)
					k = 1;
				break;
			case 15:
				if (nx <= 11)
					k = 1;
				break;
			case 16:
			case 17:
				if (nx <= 9)
					k = 1;
				break;
			case 18:
			case 19:
				if (nx <= 7)
					k = 1;
				break;
			case 20:
			case 21:
			case 22:
			case 23:
				if (nx <= 5)
					k = 1;
				break;
			default:
				if (nx <= 3)
					k = 1;
				break;
			}
#ifdef DEBUG
			status = ieee_flags("get", "exception", in, &out);
			printf("i,nx,k = %d %d %d  FLAGS = %d\n",
				i, nx, k, (status & 0x3ff));
#endif
			/* compute x**i by repeat multiplication */
			if (k == 1) {
				t = s * x;
				if (i & 1)
					z = t;
				else
					z = one;
				i >>= 1;
#ifdef DEBUG
				status =
				    ieee_flags("get", "exception", in, &out);
				printf("t = %1.20e  FLAGS = %d\n",
				    t, (status & 0x3ff));
#endif
				while (i != 0) {
					t = t * t;
					if (i & 1)
						z *= t;
					i >>= 1;
				}
#ifdef DEBUG
				status =
				    ieee_flags("get", "exception", in, &out);
				printf("t = %1.20e  FLAGS = %d\n",
				    t, (status & 0x3ff));
#endif
				return (z);
			}
		}
	}
	z = __libmcr_k_mx_pow(&x, &y, &er, &m);
#ifdef DEBUG
	printf("mx_pow: m = %d, z = %08X %08X, er = %08X %08X\n",
		m, H2(z), H2(er));
#endif
	ft = (float) rndc;
	fw = ft * ft;
	hz = HIGH_WORD(z);
	iz = hz & 0x7fffffff;
	i = (iz >> 20) + m;
	rnd = (*(int *) &fw) & 3;
	if (i > 0) {
		w = s * huge;
		if (i >= 0x7ff)
			return (w * huge);
		ncrd = __libmcr_mx_check(&z, rnd, er);
		if (ncrd == 0)
			HIGH_WORD(z) = iz + (m << 20);
	} else {	/* subnormal output, be aware of double rounding */
		w = s * tiny;
#ifdef DEBUG
		printf("subnormal output: i = %d\n", i);
#endif
		if (i < -52)
			return (w * tiny);
		HIGH_WORD(w) = ((iz >> 20) - i + 1) << 20;
		LOW_WORD(w) = (int) (tiny * tiny);
#ifdef DEBUG
		printf("w = %08X %08X %1.20e\n", H2(w), w);
#endif
		t = w + z;	/* simulate subnormal rounding by adding */
#ifdef DEBUG
		printf("t = w+z = %08X %08X %1.20e\n", H2(t), t);
		ax = t - w;
		printf("a = t-w = %08X %08X %1.20e\n", H2(ax), ax);
		ax = z - ax;
		printf("z-a = %08X %08X %1.20e\n", H2(ax), ax);
		printf("er = %08X %08X %1.20e\n", H2(er), er);
#endif
		er += z - (t - w);
#ifdef DEBUG
		printf("er + w = %08X %08X %1.20e\n", H2(er), er);
#endif
		z = t;
#ifdef DEBUG
		printf("rnd,z= %d %08X %08X %1.20e\n", rnd, H2(z), z);
		ax = z - t;
		printf("z - t= %08X %08X %1.20e\n", H2(ax), ax);
		printf("er + w+(z-t) = %08X %08X %1.20e\n", H2(er), er);
#endif
		ncrd = __libmcr_mx_check(&z, rnd, er);
		/* mask off the additional leading bit from previous adding */
		if (ncrd == 0)
			HIGH_WORD(z) &= 0x000fffff;
	}
#ifdef DEBUG
	status = ieee_flags("get", "exception", in, &out);
	printf("z = %08X %08X %1.20e %d %d\n",
		H2(z), z, ncrd, (status & 0x3ff));
	fflush(stdout);
#endif
	nw = 8;
	while (ncrd == 1) {
		z = __libmcr_mi_pow(x, y, nw, &ncrd, rnd);
#ifdef DEBUG
		status = ieee_flags("get", "exception", in, &out);
		printf("nw = %d  z = %08X %08X %1.20e %d %d\n",
			nw, H2(z), z, ncrd, (status & 0x3ff));
#endif
		nw += 2;
	}
	return (s * z);
}
