
#pragma ident "@(#)__libmcr_mm_add.c 1.3 04/02/25 SMI"

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
 * __libmcr_mm_add (int *a, int *b, int *c, int nw)
 *
 * returns c = a+b.
 * nw = no. of words for mp number
 * reliable significant bits: 24*(nw-3)
 *
 * int array:
 *	mx[0],  mx[1] ,   mx[2],   mx[3] ..... mx[nw-1]
 *	  |       |         |	   \                /
 *	sign   exponent  leading    24 bits integer
 *			  digit
 *
 * mp number structure: mi_x[0] = sign, (+1,-1)
 *			mi_x[1] = exponent,
 *			mi_x[2] = 1.0<=x<two24  (normal number)
 *			    0x7ff00000    (infinity)
 *			    0x7ff?????	(nan) (at least one ? is non-zero)
 *			    0  (zero)
 *			mp_x[3] = first significant word
 *			mp_x[4] = second significant word
 *			...
 *			mp_x[nw-1]= last significant word
 *
 * case 1: sign(x)=sign(y)
 *	1.1 exception
 *		inf/nan+inf/nan = inf/nan
 *	1.2 zero
 *		0 + b = b, a + 0 = a
 *	1.3 a>=b
 *		a[2],a[3],...,a[nw-1]
 *			  b[2],......,b[n2-1]
 *		-----------------------------
 *		c[2],c[3],....c[nw-1]
 *	1.4 a<b
 *			a[2],a[3],...,a[nw-1]
 *		b[2],........,b[nw-1]
 *		-----------------------------
 *		c[2],c[3],....c[nw-1]
 * case 2: sign(x)!=sign(y)
 *	2.1 exception
 *		inf/nan - inf/nan = nan
 *	2.2 zero
 *		0 - b = -b, a - 0 = a
 *	2.3 a>=b
 *		a[2],a[3],...,a[nw-1]
 *			  b[2],......,b[n2-1]
 *		-----------------------------
 *		  c[2],c[3],....c[nw-1]
 *	2.4 a<b
 *			a[2],a[3],...,a[nw-1]
 *		b[2],........,b[nw-1]
 *		-----------------------------
 *		c[2],c[3],....c[nw-1]
 *
 */

#include <stdlib.h>

#include "mcr.h"
#ifdef DEBUG
#include <stdio.h>
#endif

void
__libmcr_mm_add(int *a, int *b, int *c, int nw)
{
	int n, m, i, j, k, is, m24, ig1, ig2;
	int *ti, *ix, *iy;

	m24 = 0xffffff;

	if (a[0] == b[0]) {	/* true addition */
		c[0] = a[0];	/* sign of c = sign of a */
		if (a[2] >= 0x7ff00000) {	/* a is inf/nan */
			if (a[2] > b[2]) {
				c[1] = a[1];
				c[2] = a[2];
				return;
			} else {
				c[1] = b[1];
				c[2] = b[2];
				return;
			}
		} else if (b[2] >= 0x7ff00000) {	/* b is inf/nan */
			c[1] = b[1];
			c[2] = b[2];
			return;
		}
		if (a[2] == 0) {
			for (i = 1; i < nw; i++)
				c[i] = b[i];
			return;
		}		/* 0+b=b */
		if (b[2] == 0) {
			for (i = 1; i < nw; i++)
				c[i] = a[i];
			return;
		}		/* a+0=a */
		n = a[1] - b[1];
		if (n >= 0) {
			ix = a;
			iy = b;
		} else {
			ix = b;
			iy = a;
			n = -n;
		}		/* if a<b swap a and b */
		if (n >= 1 && n < nw - 1)
			k = iy[nw - n];
		else
			k = 0;
		ti = (int *) calloc(nw + 1, sizeof (int));
		for (i = nw - 1; i >= n + 2; i--) {
			ti[i + 1] = k & m24;
			k = ix[i] + iy[i - n] + (k >> 24);
		}

		for (j = i; j >= 2; j--) {
			ti[j + 1] = k & m24;
			k = ix[j] + (k >> 24);
		}
		ti[2] = k & m24;
#ifdef DEBUG
		printf("\nk, ix[2] ti[2] = %d %08X %08X\n", k, ix[2], ti[2]);
#endif
		k >>= 24;
		c[1] = ix[1] + k;  /* exponent */
		ig1 = ti[nw - k];  /* guarded word for rounding, half way up */
		if (k == 1) {	/* shift left one word */
			if (c[1] == 0x80000000) {	/* overflow */
				c[1] = 0;
				c[2] = 0x7ff00000;
			} else {
				for (i = nw - 1; i > 2; i--)
					c[i] = ti[i - 1];
				c[2] = 1;
			}
		} else
			for (i = 2; i < nw; i++)
				c[i] = ti[i];
#ifdef DEBUG
		printf("\nc[] = %d %d ", c[0], c[1]);
		for (i = 2; i < nw; i++)
			printf(" %X", c[i]);
		printf(" Guard w: %X\n", ig1);
#endif
		k = ig1 >> 23;
		i = nw - 1;
		while (k == 1 && i >= 2) {
			k += c[i];
			c[i] = k & m24;
			k >>= 24;
			i--;
		}		/* rounding */
		if (k == 1) {	/* i.e., i=1, carry one word */
			c[2] = 1;
			c[1] += 1;
		}
		free(ti);
		return;
	}
	/* end of true addition */
	else {			/* true subtraction */
		if (a[2] >= 0x7ff00000) {	/* a is inf/nan */
			if (a[2] > b[2]) {
				c[0] = a[0];
				c[2] = a[2];
				return;
			} else {
				if (b[2] == 0x7ff00000) {  /* inf - inf = nan */
					c[0] = 0;
					c[2] = b[2] + 1;
					return;
				}
				c[0] = -b[0];
				c[2] = b[2];
				return;
			}
		} else if (b[2] >= 0x7ff00000) {	/* b is inf/nan */
			c[0] = -b[0];
			c[2] = b[2];
			return;
		}
		if (a[2] == 0) { /* 0-b=-b except 0-0 = +0 */
			if (b[2] == 0)
				c[0] = 0;
			else {
				c[0] = b[0];
				for (i = 1; i < nw; i++)
					c[i] = b[i];
			} return;
		}		/* 0+b = b */
		if (b[2] == 0) {
			for (i = 0; i < nw; i++)
				c[i] = a[i];
			return;
		}		/* a-0=a */
		/* shift and subtract */
		n = a[1] - b[1];

		ti = (int *) calloc(nw + 1, sizeof (int));
		/* same exponent */
		if (n == 0) {
			/* first m that a[m]!=b[m] */
			for (m = 2; m < nw && a[m] == b[m]; m++);
			if (m == nw) {
				c[0] = 0;
				c[2] = 0;
				free(ti);
				return;
			}	/* a==b return 0 */
			if (a[m] > b[m]) {	/* if a<b swap a and b */
				ix = a;
				iy = b;
				c[1] = a[1];
				c[0] = a[0];
			} else {
				ix = b;
				iy = a;
				c[0] = b[0];
				c[1] = b[1];
			}
			k = 0x1000000;
			for (i = nw - 1; i >= m; i--) {	/* subtraction */
				k = m24 + (k >> 24) + ix[i] - iy[i];
				ti[i] = k & m24;
			}
			while (ti[m] == 0)
				m++;	/* find non-zero ti[m] */
			c[1] -= (m - 2);	/* exponent shift */
			for (i = m; i < nw; i++)
				c[i - m + 2] = ti[i];	/* normalization */
			for (j = i - m + 2; j < nw; j++)
				c[j] = 0;
			free(ti);
			return;
		}
		if (n > 0) {
			ix = a;
			iy = b;
			c[1] = a[1];
			c[0] = a[0];
		} else {
			ix = b;
			iy = a;
			c[0] = b[0];
			c[1] = b[1];
			n = -n;
		}		/* if a<b swap a and b */
		if (n >= nw) {
			for (i = 2; i < nw; i++)
				c[i] = ix[i];
			free(ti);
			return;
		}		/* ix >> iy */
		ig1 = ig2 = 0;
		is = 0;		/* 2 guard digit and sticky digit */
		/*
		 * ix[2],.......,ix[nw-1] iy[2],......................
		 * ------------------------------------ ti[2]
		 * ,ti[3],...,ti[nw-1], ig1, ig2, is
		 */
		if (n < nw - 1)
			ig1 = iy[nw - n];	/* guard digit */
		if (n >= 2) {
			if (n < nw)
				ig2 = iy[nw + 1 - n];
			if (n >= 3)
				for (i = nw + 2 - n; i < nw; i++)
					is |= iy[i];	/* sticky digit */
		}
		/* subtraction */
		if (is == 0)
			k = 0x1000000;
		else
			k = m24;
		k = k - ig2;
		ig2 = k & m24;
		k = m24 + (k >> 24) - ig1;
		ig1 = k & m24;
		for (i = nw - 1; i >= 2 + n; i--) {
			k = m24 + (k >> 24) + ix[i] - iy[i - n];
			ti[i] = k & m24;
		}
		for (j = i; j >= 2; j--) {
			k = m24 + (k >> 24) + ix[j];
			ti[j] = k & m24;
		}
		for (m = 2; ti[m] == 0; m++);	/* find non-zero ti[m] */
		if (m == 2) {
			k = ig1 >> 23;
			for (i = nw - 1; k != 0; i--) {
				k += ti[i];
				c[i] = k & m24;
				k >>= 24;
			}
			for (j = i; j >= 2; j--)
				c[j] = ti[j];
		} else if (m > 3) {
			/* happen only when n=1, do shifting (no rounding) */
			c[1] -= (m - 2);
			for (i = m; i < nw; i++)
				c[i - m + 2] = ti[m];
			c[i - m + 2] = ig1;
			for (j = i - m + 3; j < nw; j++)
				c[j] = 0;
		} else {	/* m=3, need rounding */
#ifdef DEBUG
			printf("\nti[] = ");
			for (i = 2; i < nw; i++)
				printf(" %X", ti[i]);
			printf(", g1=%X, g2=%X\n", ig1, ig2);
#endif

			k = ig1 + (ig2 >> 23);
			ig1 = k & m24;
			k >>= 24;
			for (i = nw - 1; k != 0; i--) {
				k += ti[i];
				ti[i] = k & m24;
				k >>= 24;
			}
			if (i == 1) {
				for (i = 2; i < nw; i++)
					c[i] = ti[i];
			} else {
				for (i = 2; i < nw - 1; i++)
					c[i] = ti[i + 1];
				c[nw - 1] = ig1;
				c[1] -= 1;
			}
#ifdef DEBUG
			printf("ti[] = ");
			for (i = 2; i < nw; i++)
				printf(" %X", ti[i]);
			printf(", g1=%X\n", ig1);
#endif
		}
		free(ti);
		return;
	}
}
