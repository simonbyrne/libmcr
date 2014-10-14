
#pragma ident "@(#)__libmcr_k_mx_cos.c 1.3 04/02/25 SMI"

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
 * double __libmcr_k_mx_cos(double x,double t, double *e)
 *
 * (x,t) input extended argument assuming less than ...
 * e     error of the output
 *
 * Method:
 *	1. Compute cos(x+t) by a remez approximation with precision up to
 *	   2**-79 for |x+t|<pi/4: let z = x*x,
 *		s = 1 - z/2 + (z*z) * (a1 + z* (a2 +... + z * a8)
 *
 *	Remark. starting from a4 we can use double precision arithmetic:
 *		t = z *( a4 + z * (a5 + z *(a6 + z* (a7 + z*a8))))
 *	Then we simulate higher precision for
 *		s = 1.0 - z/2 + z*z *(a1 + z*(a2 + z*(a3 + t)))
 *
 */

#include "mcr.h"

static const double
a1 = 4.16666666666666643537e-02,	/* 3FA55555 55555555 */
a1_l = 2.31293347377082407858e-18,	/* 3C455542 7F6C4C64 */
a2 = -1.38888888888888894189e-03,	/* BF56C16C 16C16C17 */
a2_l = 5.41016300000475658527e-20,	/* 3BEFEF9B 6FBDC592 */
a3 = 2.48015873015872880133e-05,	/* 3EFA01A0 1A01A016 */
a3_l = -1.30054435764941784142e-21,	/* BB98910B EEE64112 */
a4 = -2.75573192239754115370e-07,	/* BE927E4F B77897A1 */
a5 = 2.08767569835934541196e-09,	/* 3E21EED8 EFE9134E */
a6 = -1.14707445485213357474e-11,	/* BDA93974 820AFE64 */
a7 = 4.77932443029962739395e-14,	/* 3D2AE7BB 7E02157E */
a8 = -1.54972391778493833416e-16;	/* BCA65578 DB7F3581 */


double
__libmcr_k_mx_cos(x, c, e)
	double x, c, *e;
{

	double x_h, x_l, z_h, z_l, z, p, q, s, s_h, s_l, t, t_h, t_l;
	float fx, fz;
	int ix;

	ix = (HIGH_WORD(x)) & 0x7fffffff;
	if (ix < 0x3d900000) {	/* |x| < 2**(-38) */
		*e = 0.0;
		return (1.0);
	}
	fx = (float) x;
	z = x * x;
	x_h = (double) fx;
	fz = (float) z;
	x_l = c - (x_h - x);
	z_l = (x + x_h) * x_l;
	z_h = (double) fz;
	z_l += x_h * x_h - z_h;
	if (ix < 0x3ec00000) {	/* |x| < 2**(-19), use 1.0 - x*x/2 */
		fx = (float) (1.0 - 0.5 * z);
		p = (double) fx;
		q = -(0.5 * z_h - (1.0 - p));
		q -= 0.5 * z_l;
		s = p + q;
		*e = q - (s - p);
		return (s);
	}
	t = z * (a4 + z * (a5 + z * (a6 + z * (a7 + z * a8))));

	/* a3+t --> s = s_h + s_l */
	s = a3 + t;
	fx = (float) s;
	p = a3_l + t;
	s_h = (double) fx;
	t = s * z;		/* compute next t = s*z in advance */
	s_l = p - (s_h - a3);

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
	s = t - 0.5;		/* compute next s = -0.5+t in advance */
	t_l = p - (t_h - q);


	/* -0.5 + t = s = (-0.5)+(t_h+t_l) = s_h + s_l */
	fx = (float) s;
	s_l = t_l;
	s_h = (double) fx;
	t = s * z;		/* compute next t = s*z in advance */
	s_l += t_h - (s_h + 0.5);


	/* s*z = (s_h+s_l)*(z_h+z_l) = t = (t_h+t_l) */
	p = s_h * z_l + s_l * z;
	fx = (float) t;
	q = s_h * z_h;
	t_h = (double) fx;
	s = 1.0 + t;		/* compute next s = 1.0+t in advance */
	t_l = p - (t_h - q);


	/* 1.0 + t = 1.0 + t_h + t_l */
	fx = (float) s;
	s_l = t_l;
	s_h = (double) fx;
	s_l += t_h - (s_h - 1.0);
	t = s_h + s_l;
	*e = s_l - (t - s_h);


	return (t);
}
