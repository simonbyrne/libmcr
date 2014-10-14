
#pragma ident "@(#)__libmcr_k_mi_log.c 1.3 04/02/25 SMI"

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
 * void __libmcr_k_mi_log(double x, int* mc, int nw)
 *
 * Assume x is between 1/sqrt2 and sqrt2, __libmcr_k_mi_log compute
 * log(x) and return it in mi format (mc).
 *
 * Method:
 *	1. If x is really close to 1 (|r=x-1|<2**-5), then
 *		log(x) = r - 1/2 r**2 + 1/3 r**3 - ...
 *
 *	2. Otherwise, let s = (x-1)/(x+1); we have
 *		log(x)	= log(1+s)-log(1-s)
 *			= 2s + 2/3 s**3 + 2/5 s**5 + ...
 *	   Note:
 *		|s| <= (sqrt2-1)/(sqrt2+1) = 0.1715728754
 *		s*s <= 0.02943725152 ~ 2**-5.08
 *		    (i.e., every term get 5 bits; for 120 bits accuracy,
 *			   it takes ~24 terms.)
 *
 *		Let c = 2s, ss = s*s, g = 2s
 *		For i=3,5,......
 *			g = g*ss
 *			t = g/i
 *			if c[2]-t[2]>=nw stop  else
 *			c = c+g
 *
 *		How to get s in mi format? Let
 *			q0=|x-1|/(x+1)
 *		sign:		mq[0] = -1 if x<1, else mq[0] = 1
 *		exponent: 	mq[1] = -1 (initially x less than 1)
 *			q0 *= two24; xm1 = (x-1)*two24;
 *			while (q0<1.0)
 *			    { q0*=two24; xm1 *=two24;}
 *		first sig word:	mq[2] = (int)q0		... non-zero
 *
 *		Remark: mostly (int)q0 = (int) |x-1|/(x+1). It fails only
 *			when (int)|x-1|/(x+1) rounded up to an integer!
 *			In this case, 2**(mq[1]*24)*|x-1|/(x+1) = N - eps
 *			for a tiny eps. This case can be detected by
 *			examining whether the remainder is less than zero.
 *
 *		xp1_1 = x+1 chopped to 27 bits
 *		xp1_2 = (x+1)-xp1_1 (at most 27 bits)
 *		remainder = xm1 - mq[2]*(x+1)
 *			  = (xm1 - mq[2]*xp1_1) - mq[2]*xp1_2
 *		Let:
 *			r1 = two24*((x-1) - mq[2]*xp1_1) (exact);
 *			r2 = two24*(-mq[2]*xp1_2) (exact);
 *		then r = r1+r2 is the remainder.
 *		     if(r<0) {		... (int)2**(24*mq[1])*|x-1|/(x+1)
 *					...  rounded up to an int
 *			mq[2] -= 1
 *			if(mq[2]==0) {
 *			    mq[1] -= 1; mq[2] = 0xffffff;
 *			    xm1 /= two24;
 *			}
 *			r1 = two24*(xm1 - mq[2]*xp1_1);
 *			r2 = two24*(-mq[2]*xp1_2);
 *			r  = r1+r2;
 *		     }
 *
 *		Loop: for i=3,...
 *			q  = r*(1/(1+x))
 *			mq[i] = (int)q;
 *			r1 = two24*(r2-(mq[i]*xp1_1-r1));
 *			r2 = two24*(-mq[i]*xp1_2);
 *			r = r1+r2;
 *			if(r<0) 	...really rare
 *			{
 *				mq[i-1] -= 1;
 *				r1 += two24*xp1_1;
 *				r2 += two24*xp1_2;
 *			}
 *
 */

#include <stdlib.h>
#include "mcr.h"

static const double
one = 1.0,
zero = 0.0,
two24 = 16777216.0,
twom24 = 1.0 / 16777216.0;

void
__libmcr_k_mi_log(x, mc, nw)
	double x;
	int mc[], nw;
{
	double q, r, z, y, y1, y2, w, r1, r2, two24y1, two24y2;
	int j, i, *mz, *mt, *mr, *mci, nwp;
	int hx;

	nwp = nw + 1;
	mz = (int *)calloc(nwp, sizeof (int));
	mt = (int *)calloc(nwp, sizeof (int));
	mr = (int *)calloc(nwp, sizeof (int));
	mci = (int *)calloc(nwp, sizeof (int));

	hx = HIGH_WORD(x);	/* high word of x */
	if (hx >= 0x3fef0000 && hx < 0x3ff08000) {	/* |r=x-1|<2**-5 */
		/*
		 * compute log(x) by r - 1/2 *r**2 + 1/3 r**3 - ... for
		 * j=1,2,... mt = r**j, mz = r**j/j * (-1)**(j+1), mc = mc +
		 * mz.
		 */
		r = x - one;
		__libmcr_mi_dtomi(r, mr, nw);
		for (i = 0; i < nwp; i++)
			mt[i] = mci[i] = mr[i];
		j = 2;
		if (mt[1] != 0) {
			while (mci[1] - mt[1] < nwp) {
				__libmcr_mm_mul(mt, mr, mt, nwp);
				__libmcr_mm_divi(mt, j, mz, nwp);
				if ((j & 1) == 0)
					mz[0] = -1;
				__libmcr_mm_add(mci, mz, mci, nwp);
				j++;
			}
		}
	} else {		/* compute log(x) by : s=(x-1)/(x+1); log(x) */
				/* = log((1+s)/(1-s)) = 2s + 2/3 s**3 + 2/5 */
				/* s**5 + ... */
		/* compute s = |x-1|/(1+x) */
		r = x - one;
		y = x + one;
		y1 = y;
		LOW_WORD(y1) &= 0xfc000000;	/* y1 = y chopped to 27 bits */
		y2 = x - (y1 - one);	/* y2 is at most 27 bits also */
		two24y1 = two24 * y1;
		two24y2 = two24 * y2;
		if (r < zero) {
			mr[0] = -1;
			r = -r;
		} else
			mr[0] = 1;
		z = one / y;
		q = r * z;
		mr[1] = -1;
		while (q < one) {
			q *= two24;
			r *= two24;
		}
		i = (int) q;
		w = (double) i;
		mr[2] = i;
		r1 = two24 * (r - w * y1);
		r2 = two24 * (-w * y2);
		r = r1 + r2;
		if (r < zero) {	/* (int)|x-1|/(x+1) rounded up to an integer */
			mr[2] -= 1;
			if (mr[2] == 0) {
				mr[1] -= 1;
				mr[2] = 0xffffff;
				r *= twom24;
			}
			w = (double) mr[2];
			r1 = two24 * (r - w * y1);
			r2 = two24 * (-w * y2);
			r = r1 + r2;
		}
		for (i = 3; i < nwp; i++) {
			j = (int) (r * z);
			w = (double) j;
			mr[i] = j;
			r1 = two24 * (r2 - (w * y1 - r1));
			r2 = two24 * (-w * y2);
			r = r1 + r2;
			if (r < zero) {
				mr[i] = j - 1;
				r1 += two24y1;
				r2 -= two24y2;
			}
		}

		/* mt = mci = 2s initially */
		__libmcr_mm_add(mr, mr, mci, nwp);
		for (i = 0; i < nw; i++)
			mt[i] = mci[i];
		/* mr = s*s */
		__libmcr_mm_mul(mr, mr, mr, nwp);
		j = 3;
		while (mci[1] - mt[1] < nwp) {	/* if mt/mci is not too tiny */
		    __libmcr_mm_mul(mt, mr, mt, nwp);	/* mt=s*s*mt=s**j */
		    __libmcr_mm_divi(mt, j, mz, nwp);	/* mz = s**j/j */
		    __libmcr_mm_add(mci, mz, mci, nwp); /* mc=2s+...+s**j/j */
			j += 2;
		}
	}
	for (i = 0; i < nw; i++)
		mc[i] = mci[i];
	free(mz);
	free(mt);
	free(mr);
	free(mci);
}
