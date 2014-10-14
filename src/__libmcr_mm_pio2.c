
#pragma ident "@(#)__libmcr_mm_pio2.c 1.3 04/02/25 SMI"

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
 * __libmcr_mm_pio2 (int *pio2, int nw)
 * __libmcr_mm_pio4 (int *pio4, int nw)
 *
 * Return pio2 = pi/2, pio4 = pi/4, and ipio2 = 2/pi respectively
 * in nw precision.
 *
 * The algorithm that is used for computing Pi, which is due to Salamin
 * and Brent, is as follows:
 *
 * Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.
 *
 * Then from k = 1 iterate the following operations:
 *
 *     A_k = 0.5 * (A_{k-1} + B_{k-1})
 *     B_k = Sqrt (A_{k-1} * B_{k-1})
 *     D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2
 *
 * Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
 * In other words, each iteration approximately doubles the number of correct
 * digits, providing all iterations are done with the maximum precision.
 *
 */

#include <stdlib.h>
#include "mcr.h"


void
__libmcr_mm_pio2(int *p, int nw)
{
	int *ma, *mb, *md, *mt, *mw, *mp, nwp, i, k;

	if (nw < 52)
		for (i = 0; i < nw; i++)
			p[i] = __libmcr_TBL_pio2[i];
	else {
		nwp = nw + 1;
		ma = (int *) calloc(nwp, sizeof (int));
		mb = (int *) calloc(nwp, sizeof (int));
		md = (int *) calloc(nwp, sizeof (int));
		mp = (int *) calloc(nwp, sizeof (int));
		mt = (int *) calloc(nwp, sizeof (int));
		mw = (int *) calloc(nwp, sizeof (int));

		/* Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2. */
		ma[0] = mt[0] = mb[0] = 1;
		mt[1] = -1;
		mt[2] = 0x800000;
		ma[1] = mb[1] = 0;
		ma[2] = 1;
		mb[2] = 2;
		for (i = 3; i < nwp; i++)
			ma[i] = mt[i] = mb[i] = 0;
		__libmcr_mm_sqrt(mb, mb, nwp);	/* sqrt(2) */
		__libmcr_mm_sub(mb, mt, md, nwp);	/* md = sqrt(2)-1/2 */
		__libmcr_mm_scalbn(mb, nwp, -1);	/* mb = sqrt(2)/2 */

		/* loop until mp accurate to 3+24*(nw-3) bits */
		/*
		 * according to experiemnt; 1st loop get 15 bits, then
		 * 2*15+3, ...
		 */
		i = 0;

		for (k = 1; i < nw - 3; k++) {
			if (k > 3)
				i = (i << 1) + 3;
			/*
			 * A_k = 0.5 * (A_{k-1} + B_{k-1}) B_k = Sqrt
			 * (A_{k-1} * B_{k-1}) D_k = D_{k-1} - 2^k * (A_k -
			 * B_k) ^ 2
			 */
			__libmcr_mm_mul(ma, mb, mt, nwp);
			__libmcr_mm_add(ma, mb, ma, nwp);
			__libmcr_mm_scalbn(ma, nwp, -1);
			__libmcr_mm_sqrt(mt, mb, nwp);
			__libmcr_mm_sub(ma, mb, mw, nwp);
			__libmcr_mm_mul(mw, mw, mw, nwp);
			__libmcr_mm_scalbn(mw, nwp, k);
			__libmcr_mm_sub(md, mw, md, nwp);
		}
		/*
		 * P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to
		 * Pi.
		 */
		__libmcr_mm_add(ma, mb, mt, nwp);
		__libmcr_mm_mul(mt, mt, mt, nwp);
		__libmcr_mm_div(mt, md, mp, nwp);	/* mp ~ pi */
		__libmcr_mm_scalbn(mp, nwp, -1);	/* ~ pi/2 */

		for (i = 0; i < nw; i++)
			p[i] = mp[i];

		free(ma);
		free(mb);
		free(md);
		free(mp);
		free(mt);
		free(mw);
	}
}

void
__libmcr_mm_pio4(int *p, int nw)
{
	int *pio2, nwp, i;


	if (nw < 49)
		for (i = 0; i < nw; i++)
			p[i] = __libmcr_TBL_pio4[i];
	else {
		nwp = nw + 1;
		pio2 = (int *) calloc(nwp, sizeof (int));
		__libmcr_mm_pio2(pio2, nwp);
		__libmcr_mm_scalbn(pio2, nwp, -1);
		for (i = 0; i < nw; i++)
			p[i] = pio2[i];
		free(pio2);
	}
}
