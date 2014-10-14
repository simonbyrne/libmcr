
#pragma ident "@(#)__libmcr_mi_dtomi.c 1.3 04/02/25 SMI"

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

void
__libmcr_mi_dtomi(x, mx, nw)
	double x;
	int *mx, nw;
{
	int hx, lx, i, j, n;

	hx = HIGH_WORD(x);
	lx = LOW_WORD(x);
	mx[0] = ((hx >> 31) << 1) + 1;
	hx &= 0x7fffffff;
	n = (hx >> 20);
	if (n == 0) {
		if ((hx | lx) == 0) {
			for (i = 1; i < nw; i++)
				mx[i] = 0;
			return;
		}
		while (hx < 0x00100000) {
			n -= 1;
			hx = (hx << 1) + ((lx >> 31) & 1);
			lx <<= 1;
		}
		n += 1;
	}
	n -= 0x3ff;
	hx = (hx & 0xfffff) | 0x00100000;
	j = 0;
	while (n >= 24) {
		j += 1;
		n -= 24;
	}
	while (n < 0) {
		j -= 1;
		n += 24;
	}
	mx[1] = j;
	if (n <= 20) {
		if (n == 20) {
			mx[2] = hx;
			mx[3] = ((unsigned) lx) >> 8;
		} else {
			mx[2] = hx >> (20 - n);
			mx[3] = ((unsigned) (hx << (12 + n)) >> 8) |
				(((unsigned) lx) >> (28 - n));
		}
		mx[4] = (unsigned) (lx << (4 + n)) >> 8;
		if (n <= 3)
			mx[5] = (unsigned) (lx << (28 + n)) >> 8;
		else
			mx[5] = 0;
	} else {
		mx[2] = (hx << (n - 20)) | (((unsigned) lx) >> (52 - n));
		mx[3] = (unsigned) (lx << (n - 20)) >> 8;
		mx[4] = (unsigned) (lx << (4 + n)) >> 8;
		mx[5] = 0;
	}
	for (i = 6; i < nw; i++)
		mx[i] = 0;
}
