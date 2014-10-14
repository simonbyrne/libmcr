
#pragma ident "@(#)__libmcr_mm_muli.c 1.3 04/02/25 SMI"

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
 * __libmcr_mm_muli (int *a, int b, int *c, int nw)
 *
 * Assume a with nw word, b 32 bit integer, __libmcr_mm_muli returns c = a*b
 * nw = no of words for mp number
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
 * multiplication for normalize mi number :
 *		          ns significant
 *		        /      digits    \ (all converted to double)
 *	a[0],a[1],a[2], a[3],...............,a[nw-1]
 *   *	          b
 *   ----------------------------------------
 *		  b*a2, b*a3,b*a4,..........,b*a[nw-1]
 *  ===================================================
 *	t[0],t[1],t[2], t[3],t[4]...........,t[nw-1]
 *
 * Note: two extra words will guarantee 16 bits accuracy.
 *
 * Thus
 * 	ti[0]=a[0]*sign(b), 	...sign
 *	ti[1]=a[1], 	...exponent
 * t[0]=t[1]=0;
 * For i=2,...,nw-1
 *	t[i] = b*a[i]
 * If b has more than 29 bits, then break b into b1,b2, where
 * 	b1 = b>>24, b2 = b&0xffffff,
 * t[1] = b1*a[2]
 * For i=2,...,nw-1
 *	t[i] = b2*a[i-1] + b1*a[i]
 * t[nw-1] = b2*a[nw-1]
 * Then normalize t[] to ti[].
 */

#include <stdlib.h>
#include "mcr.h"

static const int
m24 = 0xffffff;

static const double
twom24 = 5.96046447753906250000e-08,
zero = 0;

void
__libmcr_mm_muli(int *a, int b, int *c, int nw)
{
	int iexp, n, m, i, j, k, *ti, Ht;
	double y1, y2, y, x, *t, da1, da2;

	ti = (int *) calloc(nw + 4, sizeof (int));
	t = (double *) calloc(nw + 4, sizeof (double));

	/* sign(c) = sign(a)*sign(b) */
	i = b >> 31;
	c[0] = (a[0] ^ ((i << 1) + 1)) + 1;
	b = (b ^ i) - i;	/* b = |b| */

	/* leading word and exception cases */
	i = a[2];
	if (i == 0 || i >= 0x7ff00000) {
		c[2] = a[2];
		free(ti);
		free(t);
		return;
	}
	/* exponent */
	iexp = a[1];

	/* compute a*b */
	if ((b & 0x1fffffff) == 0) {
		y = (double) b;
		t[1] = zero;
		k = ((nw >> 1) << 1);
		for (i = 2; i < k; i += 2) {	/* unroll 2 */
			t[i] = y * (double) a[i];
			t[i + 1] = y * (double) a[i + 1];
		}
		if (k != nw)
			t[k] = y * (double) a[k];
	} else {
		i = ((unsigned) b) >> 24;
		y1 = (double) i;
		y2 = (double) (b & 0xffffff);
		da1 = (double) a[2];
		t[1] = y1 * da1;
		for (i = 2; i < nw - 1; i++) {
			da2 = da1;
			da1 = (double) a[i + 1];
			t[i] = y2 * da2 + y1 * da1;
		}
		t[nw - 1] = y2 * da1;
	}

	/* normalize t[] to ti[] */
	t[0] = zero;
	for (i = nw - 1; i > 0; i--) {
		x = twom24 * t[i];
		Ht = HIGH_WORD(t[i]);
		m = (int) x;
		j = Ht >> 20;
		if (j == 0) {
			ti[i] = 0;
		} else {
			Ht = 0x000fffff & Ht | 0x00100000;
			n = 1075 - j;
			t[i - 1] += (double) m;
			if (n > 32)
				ti[i] = (Ht >> (n - 32)) & m24;
			else {
				k = (((unsigned) LOW_WORD(t[i])) >> n) & m24;
				if (n > 8)
					k |= ((Ht << (32 - n)) & m24);
				ti[i] = k;
			}
		}
	}
	ti[0] = m;
	ti[nw] = 0;

	/* rounding */
	if (ti[0] != 0)
		j = nw - 3;
	else if (ti[1] != 0)
		j = nw - 2;
	else
		j = nw - 1;
	k = (0x800000 + ti[j + 1]) >> 24;
	for (i = j; k != 0; i--) {
		k += ti[i];
		ti[i] = k & m24;
		k >>= 24;
	}

	/* fixing exponent and return */
	j = 0;
	if (ti[0] != 0)
		j = 2;
	else if (ti[1] != 0)
		j = 1;
	iexp += j;

	/* check for over/underflow */
	if (iexp < a[1]) {
		c[2] = 0x7ff00000;
		free(ti);
		free(t);
		return;		/* overflow */
	}
	/* return */
	c[1] = iexp;
	for (i = 2; i < nw; i++)
		c[i] = ti[i - j];
	free(ti);
	free(t);
}
