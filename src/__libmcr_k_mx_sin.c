
#pragma ident "@(#)__libmcr_k_mx_sin.c 1.3 04/02/25 SMI"

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
 * double __libmcr_k_mx_sin(double x,double t,double *e)
 *
 * (x,t) input extended argument assuming less than ...
 * e     error of the output
 *
 * Method:
 *	1. Compute sin(x+t) by a remez approximation with precision up to
 *	   2**-79 for |x+t|<pi/4: let z = x*x,
 *		s = x + x*z * (a1 + z * (a2 + ... + z * a8)
 *
 *	Remark. starting from a5 we can use double precision arithmetic:
 *		t = z * (a5 + z *(a6 + z* (a7 + z*a8)))
 *	Then we simulate higher precision for
 *		s = x*(1.0 + z *(a1 + z*(a2 + z*(a3 + z*(a4 + t)))))
 *
 */

#include "mcr.h"

static const double
a1 = -1.66666666666666657415e-01,	/* BFC55555 55555555 */
a1_l = -9.25171463293262101137e-18,	/* BC65553F 964132A9 */
a1f = -1.66666671633720397949e-01,	/* BFC55555 60000000 */
a1f_l = 4.96705373128269573703e-09,	/* 3E355555 55555603 */
a2 = 8.33333333333333321769e-03,	/* 3F811111 11111111 */
a2_l = 1.08087017178403872363e-19,	/* 3BFFE6D2 F92B1526 */
a3 = -1.98412698412698277001e-04,	/* BF2A01A0 1A01A015 */
a3_l = 1.08340943186723399741e-22,	/* 3B605F3C BB51E3F9 */
a4 = 2.75573192239740933420e-06,	/* 3EC71DE3 A556BC52 */
a4_l = 5.41023573009008118792e-23,	/* 3B5059F7 96D332F3 */
a5 = -2.50521083797723368020e-08,	/* BE5AE645 67DB1FAA */
a6 = 1.60590422489331495872e-10,	/* 3DE61245 EF0B2FD5 */
a7 = -7.64690569653402065962e-13,	/* BD6AE7B8 677CBA23 */
a8 = 2.78889974862182535120e-15;	/* 3CE91EC3 D64D9F1E */

double
__libmcr_k_mx_sin(x, c, e)
	double x, c, *e;
{
	double x_h, x_l, z_h, z_l, z, p, q, s, s_h, s_l, y, t, t_h, t_l;
	float fx, fz;
	int ix;

	ix = (HIGH_WORD(x)) & 0x7fffffff;
	if (ix < 0x3d900000) {	/* |x| < 2**(-38) */
		y = x + c;
		*e = c - (y - x);
		return (y);
	}
	fx = (float) x;
	z = x * x;
	x_h = (double) fx;
	fz = (float) z;
	x_l = c - (x_h - x);
	z_l = (x + x_h) * x_l;
	z_h = (double) fz;
	z_l += x_h * x_h - z_h;
	if (ix < 0x3ec00000) {	/* |x| < 2**(-19), use x - x*z*a1 */
		fx = (float) (1.0 + z * a1);
		p = a1f * z_h;
		s_h = (double) fx;
		q = a1f * z_l + a1f_l * (z_h + z_l);
		s_l = p - (s_h - 1.0);
		s_l += q;
		p = x_h * s_h;
		q = x_h * s_l + x_l * (s_h + s_l);
		s = p + q;
		*e = q - (s - p);
		return (s);
	}
	t = z * (a5 + z * (a6 + z * (a7 + z * a8)));

	/* a4+t --> s = s_h + s_l */
	s = a4 + t;
	fx = (float) s;
	p = a4_l + t;
	s_h = (double) fx;
	t = s * z;		/* compute next t = s*z in advance */
	s_l = p - (s_h - a4);

	/* s*z = (s_h+s_l)*(z_h+z_l) = t = (t_h+t_l) */
	p = s_h * z_l + s_l * z;
	fx = (float) t;
	q = s_h * z_h;
	t_h = (double) fx;
	s = a3 + t;		/* compute next s = a3+t in advance */
	t_l = p - (t_h - q);

	/* a3 + t = s = (a3+a3_l)+(t_h+t_l) = s_h + s_l */
	fx = (float) s;
	s_l = a3_l + t_l;
	s_h = (double) fx;
	t = s * z;		/* compute next t = s*z in advance */
	s_l += t_h - (s_h - a3);

	/* s*z = (s_h+s_l)*(z_h+z_l) = t = (t_h+t_l) */
	p = s_h * z_l + s_l * z;
	fx = (float) t;
	q = s_h * z_h;
	t_h = (double) fx;
	s = a2 + t;		/* compute next s = a2+t in advance */
	t_l = p - (t_h - q);

	/* a2 + t = s = (a2+a2_l)+(t_h+t_l) = s_h + s_l */
	fx = (float) s;
	s_l = a2_l + t_l;
	s_h = (double) fx;
	t = s * z;		/* compute next t = s*z in advance */
	s_l += t_h - (s_h - a2);

	/* s*z = (s_h+s_l)*(z_h+z_l) = t = (t_h+t_l) */
	p = s_h * z_l + s_l * z;
	fx = (float) t;
	q = s_h * z_h;
	t_h = (double) fx;
	s = a1 + t;		/* compute next s = a1+t in advance */
	t_l = p - (t_h - q);

	/* a1 + t = s = (a1+a1_l)+(t_h+t_l) = s_h + s_l */
	fx = (float) s;
	s_l = a1_l + t_l;
	s_h = (double) fx;
	t = s * z;		/* compute next t = s*z in advance */
	s_l += t_h - (s_h - a1);

	/* s*z = (s_h+s_l)*(z_h+z_l) = t = (t_h+t_l) */
	p = s_h * z_l + s_l * z;
	fx = (float) t;
	q = s_h * z_h;
	t_h = (double) fx;
	s = 1.0 + t;		/* compute next s = a1+t in advance */
	t_l = p - (t_h - q);

	/* 1.0 + t = 1.0 + t_h + t_l */
	fx = (float) s;
	s_l = t_l;
	s_h = (double) fx;
	s_l += t_h - (s_h - 1.0);

	/* s*x = (s_h+s_l)*(x_h+x_l) = t + e */
	p = s_h * x_h;
	q = s_h * x_l + s_l * x;
	t = p + q;
	s = q - (t - p);
	*e = s;
	return (t);
}
