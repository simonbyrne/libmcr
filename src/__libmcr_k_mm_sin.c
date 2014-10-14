
#pragma ident "@(#)__libmcr_k_mm_sin.c 1.3 04/02/25 SMI"

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
 * double __libmcr_k_mm_sin(mx,mc,nw)
 *
 * compute sin(mx) in mi format, where mx is assume less than pi/4
 *
 * Method:
 *	Compute sin(mr) by talyor serier
 *		mc = mx-mx**3/3!+mx**5/5!-mx**7/7!+...
 *
 *	Let mc = mx, mr = -mx*mx, mt=mx;
 *	For j=2,4,6,......
 *		mt =  mt*mr/(j*(j+1))
 *		if mc[2]-mt[2]>=nw stop  else
 *		mc = mc+mt
 *
 *	Note: M = j*(j+1) can be computed by
 *		i = 0;
 *		For j=2,4,6,...  i = i + (j<<2)-2
 *
 */

#include <stdlib.h>
#include "mcr.h"

void
__libmcr_k_mm_sin(mx, mc, nw)
	int mx[], mc[], nw;
{
	int j, i, *mt, *mr, *mci, nwp;

	nwp = nw + 1;
	mt = (int *) calloc(nwp, sizeof (int));
	mr = (int *) calloc(nwp, sizeof (int));
	mci = (int *) calloc(nwp, sizeof (int));

	/* Taylor series */
	for (i = 0; i < nw; i++)
		mt[i] = mci[i] = mx[i];
	for (i = nw; i < nwp; i++)
		mt[i] = mci[i] = 0;
	__libmcr_mm_mul(mt, mt, mr, nwp);
	mr[0] = -1;		/* mr = -mx*mx */
	if (mr[2] != 0) {
		j = 2;
		i = 0;		/* i = j*(j+1) */
		while (mc[1] - mt[1] < nwp - 1) {  /* mt = mt*r/j = r**j/j! */
			__libmcr_mm_mul(mt, mr, mt, nwp);
			i += (j << 2) - 2;
			__libmcr_mm_divi(mt, i, mt, nwp);  /* mt/(j*(j+1)) */
			__libmcr_mm_add(mci, mt, mci, nwp);
			j += 2;
		}
	}
	for (i = 0; i < nw; i++)
		mc[i] = mci[i];

	free(mt);
	free(mr);
	free(mci);
}
