
#pragma ident "@(#)__libmcr_mi_itomi.c 1.3 04/02/25 SMI"

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
 * mi_itomi (int i, int *mx, int nw)
 *
 * Convert i to mx number.
 * nw = no of words for mp number
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
 */

#include <stdlib.h>
#include "mcr.h"

static const int
m24 = 0xffffff;

void
__libmcr_mi_itomi(int i, int *mx, int nw)
{
	int k;

	/* sign(mx) = sign(i) */
	k = i >> 31;
	mx[0] = k | 1;
	i = (i ^ k) - k;	/* abs(i) */

	/* exponent and leading word */
	if ((i & 0xff000000) == 0) {
		mx[1] = 0;
		mx[2] = i;
		mx[3] = 0;
	} else {
		mx[1] = 1;
		mx[2] = ((unsigned) i) >> 24;
		mx[3] = i & m24;
	}

	for (i = 4; i < nw; i++)
		mx[i] = 0;
}
