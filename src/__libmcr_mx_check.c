
#pragma ident "@(#)__libmcr_mx_check.c 1.7 04/02/25 SMI"

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

#ifdef DEBUG
#include <stdio.h>
#endif

int
__libmcr_mx_check(z, rnd, er)
	double *z, er; int rnd;
/*
 * er is the error correction of z, with at least 20 bit accuracy. Assume
 * z=(iz,lz),er=(hr,0) is normal. Check whether er is close to 1/2 ulp (if
 * rnd = 0 -- Round-to-Nearest) or 0 or an ulp (if rnd !=0 --
 * Round-to-zero,inf,-inf)
 */
{
	int ncrd, iz, hz, ie, lz, hr, j, k;
	double w;
	hz = HIGH_WORD(*z);
	lz = LOW_WORD(*z);
	hr = HIGH_WORD(er);

#ifdef DEBUG
	printf("__libmcr_mx_check: hz,lz,hr,rnd = %08X %08X %08X %d\n",
		hz, lz, hr, rnd);
#endif

	/* no subnormal result */
	ie = hr & 0x7fffffff;
	iz = hz & 0x7ff00000;

	/* determine correctly roundingness by asuming 20 bit extra accuracy */
	ncrd = 0;
	if (rnd == 0) {		/* rounded to nearest */
		iz -= 0x03500000;	/* half ulp of z */
		if ((lz | (hz & 0x000fffff)) == 0 && ((hz ^ hr) < 0))
			iz -= 0x00100000;
		j = 4 + ((ie & 0x00000100) >> 6);
		k = (ie + j) & 0xfffffff8;
		if (k == iz)
			ncrd = 1;
		else if (k > iz) {	/* error > 0.5ulp, need adjustment */
			if ((hz ^ hr) >= 0) {
				k = LOW_WORD(*z) + 1;
				if (k == 0)
					HIGH_WORD(*z) += 1;
				LOW_WORD(*z) = k;
			} else {
				k = LOW_WORD(*z);
				if (k == 0)
					HIGH_WORD(*z) -= 1;
				LOW_WORD(*z) = k - 1;
			}
		}
	} else {		/* rounded toward zero, inf, or -inf */
		/* final adjustment of direct rounding */
		w   = *z + er;
		er += *z - w;
		*z  = w;
		iz -= 0x03400000;	/* 1 ulp of z */
		if ((lz | (hz & 0x000fffff)) == 0 && ((hz ^ hr) < 0))
			iz -= 0x00100000;
		j = 2 + ((ie & 0x00000100) >> 7);
		k = (ie + j) & 0xfffffffc;
		j = (iz - ie) >> 20;
		if (k == iz || k == (iz + 0x00100000) || j >= 19)
			ncrd = 1;	/* tiny error  <= 2**-20 ulp */
	}
#ifdef DEBUG
	printf("   .....: iz,ie,rnd = %08X %08X %d\n", iz, ie, rnd);
	printf("ncrd = %d\n", ncrd);
#endif
	return (ncrd);
}
