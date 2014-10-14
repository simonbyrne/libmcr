
#pragma ident "@(#)__libmcr_mm_scalbn.c 1.3 04/02/25 SMI"

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

#include "mcr.h"

/*
 * void __libmcr_mm_scalbn(mx,nw,n)
 * compute mi*2**n by scaling the exponent
 */

void
__libmcr_mm_scalbn(mx, nw, n)
	int mx[], nw, n;
{
	int i, j, jk, k, m;

	i = mx[2] & 0x7fffffff;

	/* 0, inf or nan * 2**n is still 0, inf or nan */
	if (i == 0 || i >= 0x7ff00000)
		return;

	/* normalize n to [0,24) */
	while (n >= 24) {
		n -= 24;
		mx[1] += 1;
	}
	while (n < 0) {
		n += 24;
		mx[1] -= 1;
	}

	/* compute k = logb(mx[2]) */
	k = 1;
	while ((i >>= 1) != 0)
		k += 1;

	/*
	 * if k+n>24, shift right (24-n). Note that since k+n>=24, k>=24-n
	 */

	m = 24 - n;
	if (n == 0)
		return;
	if (k + n > 24) {
		mx[1] += 1;
		j = mx[2];
		mx[2] = j >> m;
		k = 32 - m;
		for (i = 3; i < nw; i++) {
			jk = mx[i];
			mx[i] = (jk >> m) | (((unsigned) (j << k)) >> 8);
			j = jk;
		}
	} else {

		/* if k+n<24, shift left n. */
		j = mx[2];
		jk = mx[3];
		mx[2] = (j << n) | ((jk >> m));
		j = jk;
		for (i = 3; i < nw - 1; i++) {
			jk = mx[i + 1];
			mx[i] = ((j << n) & 0xffffff) | ((jk >> m));
			j = jk;
		}
		mx[nw - 1] = (jk << n) & 0xffffff;
	}
}
