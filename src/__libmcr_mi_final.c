
#pragma ident "@(#)__libmcr_mi_final.c 1.4 04/02/25 SMI"

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
 * double __libmcr_mi_final(int mc[],int nw,int *ncrd,int rnd)
 *
 * Convert mi format mc into a double with respect to the
 * rounding dirction rnd with *ncrd=1 to indicate the result is possibly
 * not correctly rounded (hence recomputation is needed). Here nw is the
 * number of integer words in the multi-precision array.
 *
 * According to our analysis, all basic routines in mm*.c produce results
 * with 24*(nw-3) significant bits. It is thus safe to allow UNC(=10) bits of
 * uncertainty for routine in __libmcr_mi_*.c.
 *
 * To determine whether the error is close to an ulp or 0.5 ulp,
 * four values will be given by the subroutine __libmcr_mi_mitod:
 *
 *    |---------------------------------|
 *    |	1.xxx ....xx   r   ss........s  |
 *    |---------------------------------|
 *
 *     bit        bit bit bit       bit
 *      0         53  54  55 -> 24(nw-3)-UNC
 *                   (mrd) (mst,ms0,ms1)
 * mtail[0] = mrd --- round bit (the 1st-bit after 53rd bit)
 * mtail[1] = mst --- sticky bit (if all s are 0, then 0; otherwise 1)
 * mtail[2] = ms0 --- =1 if all s are 0 (= ~mst)
 * mtail[3] = ms1 --- =1 if all s are 1
 *
 * Thus, assuming the precision is good up to 24*(nw-2)-UNC, we can determine
 * whether the result is correctly rounded (or chopped) as follow.
 * Check whether the answer is near half ulp
 *		*ncrd = 0;
 *		if(mrd==1&&ms0==1) *ncrd = 1;
 *		if(mrd==0&&ms1==1) *ncrd = 1;
 * Check whether the answer is near an ulp
 *		*ncrd = 0
 *		if(mrd==1&&ms1==1) *ncrd = 1;
 *		if(mrd==0&&ms0==1) *ncrd = 1;
 */

#include "mcr.h"
#ifdef DEBUG
#include <stdio.h>
#endif

#define	MGD mtail[0]		/* guard bit */
#define	MST mtail[1]		/* sticky bit */
#define	MS0 mtail[2]		/* =1 if all bits after the guard bit up to */
				/* the last UNC are 0 */
#define	MS1 mtail[3]		/* =1 if all bits after the guard bit up to */
				/* the last UNC are 1 */

double
__libmcr_mi_final(mc, nw, ncrd, rnd)
	int mc[], nw, *ncrd, rnd;
{
	/*
	 * rnd: rounding direction: 0 --  nearest, 1 - chopped, 2 - up, 3 -
	 * down
	 */
	double r;
	int xsb, mtail[4];
	int hx, lx;

	*ncrd = 0;

	/* convert to double */
	r = __libmcr_mi_mitod(mc, nw, mtail);

	hx = HIGH_WORD(r);
	lx = LOW_WORD(r);
	xsb = hx & 0x80000000;
#ifdef DEBUG
	printf("MC = % d %d %X %X %X .... %X %X\n",
	    mc[0], mc[1], mc[2], mc[3], mc[4], mc[5], mc[nw - 2], mc[nw - 1]);
	printf("to double: %08X %08X %d %d\n",
		HIGH_WORD(r), LOW_WORD(r), MGD, MST);
	printf("Tail all 0? %d\n", MS0);
	printf("Tail all 1? %d\n", MS1);
#endif

	/* rounding according to rnd */
	switch (rnd & 3) {
	case 0:		/* nearest */
		if ((MGD == 1) && (MST != 0 || (lx & 1) != 0)) {
			lx += 1;
			if (lx == 0)
				hx += 1;
		}
		if ((MGD == 1 && MS0 == 1) || (MGD == 0 && MS1 == 1))
			*ncrd = 1;
		break;
	case 1:		/* chopped */
		if ((MGD == 1 && MS1 == 1) || (MGD == 0 && MS0 == 1))
			*ncrd = 1;
		break;
	case 2:		/* to +inf */
		if ((MGD == 1 && MS1 == 1) || (MGD == 0 && MS0 == 1))
			*ncrd = 1;
		if (xsb == 0 && (MGD | MST) != 0) {
			lx += 1;
			if (lx == 0)
				hx += 1;
		}
		break;
	case 3:		/* to -inf */
		if ((MGD == 1 && MS1 == 1) || (MGD == 0 && MS0 == 1))
			*ncrd = 1;
		if (xsb != 0 && (MGD | MST) != 0) {
			lx += 1;
			if (lx == 0)
				hx += 1;
		}
		break;
	}

	/* return */
	HIGH_WORD(r) = hx;
	LOW_WORD(r) = lx;
	return (r);
}
