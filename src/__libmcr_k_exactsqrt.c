
#pragma ident "@(#)__libmcr_k_exactsqrt.c 1.2 04/02/25 SMI"

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

#include <math.h>	/* sqrt */
#include "mcr.h"

/*
 * __libmcr_k_exactsqrt(double z, int *inx) computes sqrt(z)
 * and indicates whether the result is exact (*inx = 0 -- exact,
 * *inx = 1, inexact).
 * Method 1:
 *	Let w = sqrt(z), nz = ilogb(z) and nw = ilogb(w). If mw
 *	denote the length of non-zero mantissa of w, then w is
 *	the exact sqrt of z if and only if (mw <=26 or
 *	(mw=27 and nw+nw=nz)) and w*w=z.
 * Method 2:
 *	First, clear the accrued inexact flag. Next Compute
 *	z = sqrt(z). Finally check the accrued inexact flag
 *	to see if the sqrt is exact or not, and set the value
 *	of *inx accordingly.
 */

#if !defined(__sparc)	/* Method 1 */
double
__libmcr_k_exactsqrt(double z, int *inx)
{
	double w;
	int i, j, lw, hw;
	*inx = 0;
	w = sqrt(z);
	i = (HIGH_WORD(z) >> 20) - 1023;
	lw = LOW_WORD(w);
	hw = HIGH_WORD(w);
	if ((lw & 0x3ffffff) == 0) {
		if ((lw & 0x4000000) != 0) {
			j = (hw >> 20) - 1023;
			if (j + j != i)
				*inx = 1;
			else if (w * w != z)
				*inx = 1;
		}
	} else
		*inx = 1;
	return (w);
}

#else	/* Method 2 */
#include <ieeefp.h>
extern fp_except fpsetsticky(fp_except);	/* change logged exceptions */
extern fp_except fpgetsticky(void);		/* return logged exceptions */

double
__libmcr_k_exactsqrt(double z, int *inx)
{
	fp_except ex;
	*inx = 0;
	fpsetsticky(FP_CLEAR);
	z = sqrt(z);
	ex = fpgetsticky();
	if (ex == FP_SET)
		*inx = 1;
	return (z);
}
#endif
