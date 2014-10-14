
#pragma ident "@(#)__libmcr_mm_mul.c 1.3 04/02/25 SMI"

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
 * __libmcr_mm_mul (int *a, int *b, int *c, int nw)
 *
 * Assume a,b with nw word, __libmcr_mm_mul returns c = a*b
 * nw = no of words for mp number
 * Accuracy: no of the reliable digits: 24*(nw-3)
 *
 * int precision array:
 *	mx[0],  mx[1] ,   mx[2],   mx[3] ..... mx[nw-1]
 *	  |       |         |	   \                /
 *	sign    binary,  leading    24 bits integer
 *	       exponent   digit
 *
 * mi number structure: mi_x[0] = sign, (+1,-1)
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
 *	a[0],a[1],a[2], a[3],.............,a[nw-1]
 *   *	b[0],b[1],b[2], b[3],.............,b[nw-1]
 *   --------------------------------------
 *		  b2a2, b2a3,b2a4,........,b2a[nw-1]
 *		        b3a2,b3a3,........,b3a[nw-2],b3a[nw-1]
 *		             b4a2,........,b4a[nw-3],b4a[nw-2],b4a[nw-1]
 *				  ..
 *		+			   b[nw-1]a2,b[nw-1]*a3,b[nw-1]*a4,...
 *  ============================================================================
 *	t[0],t[1],t[2], t[3],t[4].........,t[nw-1],  t[nw],  t[nw+1],   t[nw+2]
 *
 * Note 1.  If nw < 261,
 *	the two extra words in the final sum will guarantee 16 bits
 *      accuracy in t[nw]. Three extra words will allow nw < 2**32+5 for the
 *      same accuracy.
 *	Proof: the t[nw+2] will be log2(nw-4)+48 wide and hence will affect
 *	       t[nw] by log2(nw-4) bits. Now log2(nw-4)<8 => nw < 261.
 *	       Similarly, t[nw+3] will be log2(nw-5)+48 wide. It will affect
 *	       t[nw] by log2(nw-5)-24. If log2(nw-5)-24 < 8, we will have
 *	       log2(nw-5)<32 ==> nw-5 < 2**32  ==> nw < 2**32+5. QED.
 *	(We will compute 3 extra words.)
 * Note 2.Each product has 48 bits and each IEEE double has 53 significant bits.
 *	Thus each IEEE double can hold the sum of 32 products exactly.
 *	Special care is needed when summing more than 32 products.
 *
 * Thus (assume k=0 and t[i]=0 for i=0,...,nw+2)
 * 	t[0]=a[0]*b[0], 	...sign
 *	t[1]=a[1]+b[1], 	...exponent
 *
 * For i=2,...,nw-1
 *	k += 1
 *	t[j] += b[i]*a[2+(j-i)] for j=i,...,min{nw+2,nw+(i-3)}
 *	if k>=30 ... more than 30 term added, need to normalize t[i]
 *		     to avoid inexact arithmetic
 *	    for(jj=nw+2;jj>i;j--) {
 *		    x = (double)((int)(twom48*t[jj]));
 *		    t[jj]   -= x*two48;
 *		    t[jj-2] += x;
 *	    }
 *	k = 0;		... reset k to zero
 *
 */
#include <stdlib.h>
#include "mcr.h"

static const double
twom24 = 5.96046447753906250000e-08,	/* 2**(-24) */
twom48 = 3.55271367880050092936e-15,	/* 2**(-48) */
two48 = 281474976710656.0,	/* 2**(+48) */
two24 = 16777216.0,		/* 2**(+24) */
zero = 0;

void
__libmcr_mm_mul(int *a, int *b, int *c, int nw)
{
	int iexp, n, m, i, j, jj, *ti;
	double z, y, x, *ta, *t;

	ta = (double *) calloc(nw + 4, sizeof (double));
	t = (double *) calloc(nw + 4, sizeof (double));
	ti = (int *) calloc(nw + 4, sizeof (int));

	/* sign(c) = sign(a)*sign(b) */
	c[0] = (a[0] ^ b[0]) + 1;

	/* leading word and exception cases */
	i = a[2];
	j = b[2];
	if (i == 0 || i >= 0x7ff00000 || j == 0 || j >= 0x7ff00000) {
		if (i == 0) {
			if (j >= 0x7ff00000) {
				if (j > 0x7ff00000)
					c[2] = b[2];
				else
					c[2] = 0x7fffffff;
			} else
				c[2] = 0;
		} else if (j == 0) {
			if (i >= 0x7ff00000) {
				if (i > 0x7ff00000)
					c[2] = b[2];
				else
					c[2] = 0x7fffffff;
			} else
				c[2] = 0;
		} else {
			if (i > j)
				c[2] = a[2];
			else
				c[2] = b[2];
		}
		for (i = 3; i < nw; i++)
			c[i] = 0;
		c[1] = 0;
		free(t);
		free(ti);
		free(ta);
		return;
	}
	/* exponent */
	iexp = a[1] + b[1];

	/* compute the prouct of a*b */

	/* initialize t[i] = b[2]*a[i]+b[3]*a[i-1]+b[4]*a[i-2] */
	x = (double) b[2];
	y = (double) b[3];
	z = (double) b[4];
	ta[2] = (double) a[2];
	ta[3] = (double) a[3];
	t[2] = x * ta[2];
	t[3] = x * ta[3] + y * ta[2];

	for (i = 4; i < nw; i++) {
		ta[i] = (double) a[i];
		t[i] = x * ta[i] + y * ta[i - 1] + z * ta[i - 2];
	}
	t[nw] = y * ta[nw - 1] + z * ta[nw - 2];
	t[nw + 1] = z * ta[nw - 1];
	t[nw + 2] = zero;

	n = 3;
	for (j = 5; j < nw; j++) {
		if (b[j] != 0) {
			n += 1;
			y = (double) b[j];
			m = 2 - j;
			for (i = j; i <= nw + 2; i += 1)
				t[i] += y * ta[m + i];
		}
		if (n > 30) {
			/*
			 * more than 30 term have been added, need to
			 * normalize t[i] to avoid inexact arithmetic
			 */
			for (jj = nw + 2; jj > j; jj--) {
				x = (double) ((int) (twom48 * t[jj]));
				t[jj] -= x * two48;
				t[jj - 2] += x;
			}
			n = 0;	/* reset n to zero */
		}
	}

	/* normalize t[] to ti[] */
	t[0] = t[1] = zero;
	for (i = nw + 2; i > 0; i--) {
		x = (double) ((int) (twom24 * t[i]));
		ti[i] = (int) (t[i] - x * two24);
		t[i - 1] += x;
	}
	ti[0] = (int) t[0];

	/* fixing exponent */
	j = 0;
	if (ti[0] != 0)
		j = 2;
	else if (ti[1] != 0)
		j = 1;
	iexp += j;

	/* check for over/underflow */
	if (a[1] >= b[1]) {
		if (a[1] >= 0) {
			if (b[1] >= 0 && a[1] > iexp) {	/* overflow */
				c[2] = 0x7ff00000;
				free(t);
				free(ti);
				free(ta);
				return;
			}
		} else if (iexp > b[1]) {	/* underflow */
			c[2] = 0x0;
			free(t);
			free(ti);
			free(ta);
			return;
		}
	} else {
		if (b[1] >= 0) {
			if (b[1] > iexp && a[1] >= 0) {	/* overflow */
				c[2] = 0x7ff00000;
				free(t);
				free(ti);
				free(ta);
				return;
			}
		} else if (iexp > a[1]) {	/* underflow */
			c[2] = 0x0;
			free(t);
			free(ti);
			free(ta);
			return;
		}
	}

	/* return */
	c[1] = iexp;
	for (i = 2; i < nw; i++)
		c[i] = ti[i - j];
	free(t);
	free(ti);
	free(ta);
}
