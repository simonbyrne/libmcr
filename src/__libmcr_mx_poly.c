
#pragma ident "@(#)__libmcr_mx_poly.c 1.3 04/02/25 SMI"

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
 * (void) __libmcr_mx_poly (double *z, double *a, double *e, int n)
 * return
 *	e = a  + z*(a  + z*(a  + ... z*(a  + e)...))
 *	     0       2       4           2n
 * Note:
 * 1.	e and coefficient ai are represented by two double numbers.
 *	For e, the first one contain the leading 24 bits rounded, and the
 *	second one contain the remaining 53 bits (total 77 bits accuracy).
 *	For ai, the first one contian the leading 53 bits rounded, and the
 *	second is the remaining 53 bits (total 106 bits accuracy).
 * 2.	z is an array of three doubles.
 * 	z[0] :	the rounded value of Z (the intended value of z)
 * 	z[1] :	the leading 24 bits of Z rounded
 * 	z[2] :	the remaining 53 bits of Z
 *		Note that z[0] = z[1]+z[2] rounded.
 *
 */

#include "mcr.h"

void
__libmcr_mx_poly(z, a, e, n)
	double z[], a[], e[];
	int n;
{
	double r, s, t, p_h, p_l, z_h, z_l, p;
	int i;

	n = n + n;
	p = e[0] + a[n];
	p_l = a[n + 1];
	p_h = (double) ((float) p);
	p = a[n - 2] + z[0] * p;
	z_h = z[1];
	z_l = z[2];
	p_l += e[0] - (p_h - a[n]);

	for (i = n - 2; i >= 2; i -= 2) {

		/* compute p = ai + z * p */
		t = z_h * p_h;
		s = z[0] * p_l + p_h * z_l;
		p_h = (double) ((float) p);
		s += a[i + 1];
		r = t - (p_h - a[i]);
		p = a[i - 2] + z[0] * p;
		p_l = r + s;
	}
	e[0] = (double) ((float) p);
	t = z_h * p_h;
	s = z[0] * p_l + p_h * z_l;
	r = t - (e[0] - a[0]);
	e[1] = r + s;
}
