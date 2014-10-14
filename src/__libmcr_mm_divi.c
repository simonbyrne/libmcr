
#pragma ident "@(#)__libmcr_mm_divi.c 1.3 04/02/25 SMI"

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
 * __libmcr_mm_divi (int *a, int b, int *c, int nw)
 *
 * Assume a with nw word, b is an integer,  __libmcr_mm_divi
 * returns c = a[]/b chopped
 * nw = no of words for mp number
 * Reliable significant bits: 24*(nw-3)
 *
 * int precision array:
 *	mx[0],  mx[1] ,   mx[2],   mx[3] ..... mx[nw-1]
 *	  |       |         |	   \                /
 *	sign    binary,  leading    24 bits integer
 *	       exponent   digit
 *
 * mi number structure: mi_x[0] = sign,
 *			mi_x[1] = exponent (binary),
 *			mi_x[2] = 1.0<=x<two24  (normal number)
 *				  0x7ff00000	(infinity)
 *				  >0x7ff00000	(nan)
 *				  0 		(zero)
 *			mi_x[3] = first significant word
 *			mi_x[4] = second significant word
 *			...
 *			mi_x[nw-1]= last significant word
 *
 * division for normalize mi number :
 *	    t[2],t[3],.......
 *	   ______________________________________________
 *	b / a[2],a[3],a[4], a[5],...............,a[nw-1]
 *
 * Note: to get a correctly rounded result, we need to calculate up to t[nw+2].
 *	For if a[2]=1 and b is 32 bits, then t[2]=t[3]=0.
 *	So two guard words and one sticky word are needed.
 *
 * Thus (assume t[i]=0 for i=0,...,nw+2)
 * 	t[0]=a[0]*sign(b), 	...sign
 *	t[1]=a[1], 		...exponent
 *	db = (double)|b|.
 * If b is less than or equal to 29 bits, then t[i]*b is exact
 * (because t[i] is at most 24 bits). Thus
 *   Let r=(double)a[2].
 *   For i=2,...,nw+2 {
 *	t[i]   = (int)(r/db); 		... should be computed in chopped mode.
 *	tt = r-db*(double)t[i];
 *	r  = 2**24*tt;
 *   		... if r/db is not computed in chopped mode, then tt may not be
 *		... positive.
 *	if (tt<0) {
 *	    r += 2**24*db; t[i] -= 1;
 *	}
 *	if(i<nw-1) r += (double)a[i+1];
 *   }
 *
 * If b is more than 29 bits.
 *   Let b1 be the lower 29 bits of b and b0 = b-b1.
 *   Let r = (double)a[2]*2**24, t[2] = 0
 *   For i=3,...,nw+2 {
 *	t[i]   = (int)(r/db);
 *	dt = (double)t[i];
 *	tt = ((r-(b0*dt))-b1*dt);
 *	if(i<nw&&a[i]!=0) tt += (double)a[i];
 *	r  = 2**24*tt;
 *	if (tt<0) {
 *	    r += 2**24*db; t[i] -= 1;
 *	}
 *	if (tt>db) {
 *	    r -= 2**24*db; t[i] += 1;
 *	}
 *   }
 *
 */

#include <stdlib.h>
#include "mcr.h"

#ifdef DEBUG
#include <stdio.h>
#endif

static const double
zero = 0.0,
two24 = 16777216.0;

void
__libmcr_mm_divi(int *a, int b, int *c, int nw)
{
	int iexp, m, i, j, k, *ti;
	double r, y, y0, y1, x, w, tt, xj;

	/* sign(c) = sign(a)*sign(b) */
	i = b >> 31;
	c[0] = (a[0] ^ i) | 1;
	/* b = |b|, note that b>=0 except when b=0x80000000 */
	b = (b ^ i) - i;

	/* leading word and exception cases */
	i = a[2];
	if (i == 0 || i >= 0x7ff00000) {
		c[2] = a[2];
		c[1] = 0;
		for (j = 2; j < nw; j++)
			c[j] = 0;
		if ((i | b) == 0)
			c[2] = 0x7fffffff;
		return;
	}
	if (b == 0) {
		c[2] = 0x7fffffff;
		c[1] = 0;
		for (j = 2; j < nw; j++)
			c[j] = 0;
		return;
	}
	ti = (int *) calloc(nw + 4, sizeof (int));

	/* exponent */
	iexp = a[1];

	/* compute the quotient of a[]/b */
	y = (double) b;
	if (b < 0)
		y = -y;
	w = y * two24;
	k = b & 0xe0000000;
	if (k == 0) {
		r = (double) a[2];
		for (i = 2; i <= nw + 2; i++) {
			x = r / y;
			j = (int) x;
			tt = r - y * (double) j;
			r = two24 * tt;
			if (tt < zero) {
				r += w;
				j -= 1;
			}
			if (i < nw - 1)
				r += (double) a[i + 1];
			ti[i] = j;
		}
	} else {
		y0 = (double) k;
		if (b < 0)
			y0 = -y0;
		i = b & 0x1fffffff;
		y1 = (double) i;
		ti[2] = 0;
		r = two24 * (double) a[2];
		for (i = 3; i <= nw + 2; i++) {
			x = r / y;
			j = (int) x;
			xj = (double) j;
			tt = (r - y0 * xj) - y1 * xj;
			if (i < nw - 1 && a[i] != 0)
				tt += (double) a[i];
			r = two24 * tt;
			if (tt < zero) {
				r += w;
				j -= 1;
			} else if (tt > y) {
				r -= w;
				j += 1;
			}
			ti[i] = j;
		}
	}
	/* determine the sticky word */
#ifdef DEBUG
	for (i = 0; i <= nw + 2; i++)
		printf("__libmcr_mm_divi q[%2d] = %08X\n", i, ti[i]);
#endif
	if (ti[2] == 0) {
		if (ti[3] == 0) {
			m = nw + 1;
			iexp -= 2;
		} else {
			m = nw;
			iexp -= 1;
		}
	} else
		m = nw - 1;

	/* shifting (no rounding, chopped) */
	for (i = nw - 1; i >= 2; m--, i--)
		c[i] = ti[m];

	free(ti);

	/* check for over/underflow */
	if (a[1] < iexp) {
		c[2] = 0x0;
		return;		/* underflow */
	}
	c[1] = iexp;
}
