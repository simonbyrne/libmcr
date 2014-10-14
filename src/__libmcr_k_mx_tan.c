
#pragma ident "@(#)__libmcr_k_mx_tan.c 1.3 04/02/25 SMI"

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
 * double __libmcr_k_mx_tan(double x,double t,double *e)
 *
 * (x,t) input extended argument assuming less than ...
 * e     error of the output
 *
 * Method:
 *	1. Compute tan(x+t) by a remez approximation with precision up to
 *		2**-78 for |x+t|<pi/4: let z = x*x,
 *		s = x + x*z * (a1 + z * (a2 + ... + z * a20)
 *
 *	Remark. starting from a13 we can use double precision arithmetic:
 *		t = z * (a13 + z *(a14 + ... z* (a19 + z*a20)))
 *	Then we simulate higher precision for
 *		s = x*(1.0 + z *(a1 + z*(a2 + z*(a3 + ...  z*(a12 + t)))))
 *
 */

static double a[] = {
	1.0,
	0.0,
	3.33333333333333314830e-01,	/* a1   3FD55555 55555555 */
	1.85020658985619392666e-17,	/* a1_  3C7554D8 92E54E4C */
	1.33333333333333331483e-01,	/* a2   3FC11111 11111111 */
	2.32141804857483275473e-18,	/* a2_  3C45694B 1BC80441 */
	5.39682539682539222370e-02,	/* a3   3FABA1BA 1BA1BA15 */
	-1.36842763754436593729e-18,	/* a3_  BC393E37 80D86037 */
	2.18694885361576385474e-02,	/* a4   3F9664F4 882C13B8 */
	4.30260692360529407207e-19,	/* a4_  3C1FBF65 0A675C13 */
	8.86323552982673339151e-03,	/* a5   3F8226E3 55E6184F */
	8.34212077893300381996e-19,	/* a5_  3C2EC6E9 0A3D6295 */
	3.59212803811879630156e-03,	/* a6   3F6D6D3D 0E4BE5D8 */
	-1.73249639421640312385e-19,	/* a6_  BC09912F B2DD1555 */
	1.45583436480701622698e-03,	/* a7   3F57DA36 3F0E282A */
	-1.19222330001090716403e-20,	/* a7_  BBCC268C 32D5314B */
	5.90027674793155943941e-04,	/* a8   3F435582 C89285E6 */
	-2.61478629881606844265e-20,	/* a8_  BBDEDEB4 E88E8274 */
	2.39127266851269355786e-04,	/* a9   3F2F57C7 94D7B00F */
	-2.96455154230951050154e-22,	/* a9_  BB766646 D2339D33 */
	9.69265612770624793507e-05,	/* a10  3F1968A1 A4A86BA8 */
	-3.02648068762705435590e-21,	/* a10_ BBAC9594 CF87FC8F */
	3.92257978345854642333e-05,	/* a11  3F0490CC 269F9323 */
	2.85582089367268506137e-21,	/* a11_ 3BAAF8F3 83C0459D */
	1.61118410817620263042e-05,	/* a12  3EF0E4FD 49D4137C */
	1.15036173557066262674e-21,	/* a12_ 3B95BACE DA24D8E5 */
	5.89613724877126392884e-06,	/* a13  3ED8BAED DDB1F75D */
	-1.42200395896102832139e-22,	/* a13_ BB657D1C 826BB7CB */
	3.86557292923440222823e-06,	/* a14  3ED036A0 988EB1B2 */
	3.85668077595137942625e-22,	/* a14_ 3B7D23E7 D20D666A */
	-1.12346990848353174707e-06,	/* a15  BEB2D944 3A53559A */
	-9.88494678968876892106e-23,	/* a15_ BB5DE01B 6D3F365A */
	3.33546612615727943157e-06,	/* a16  3ECBFADB E516F360 */
	7.20602812664093892358e-23,	/* a16_ 3B55C764 EB332DAB */
	-2.69304350848159810452e-06,	/* a17  BEC69744 536A5273 */
	1.60298770520888692132e-22,	/* a17_ 3B683942 21555361 */
	2.06618667835868115171e-06,	/* a18  3EC1551A 2395DF50 */
	2.08421399434678964974e-22,	/* a18_ 3B6F7EE9 903003B7 */
	-8.62996260922438491137e-07,	/* a19  BEACF514 D8E398FD */
	5.23514677147119747097e-23,	/* a19_ 3B4FA4FF 4BC564AA */
	2.14535976383859239829e-07,	/* a20  3E8CCB66 67A66395 */
	1.22630029761662771020e-23,	/* a20_ 3B2DA66E 62AAB3EB */
};

#include "mcr.h"

double
__libmcr_k_mx_tan(x, c, err)
	double x, c, *err;
{

	double x_h, x_l, z_h, z_l, r, z, s, y, t, ee[2], zz[3];
	int ix, i;

	ix = (HIGH_WORD(x)) & 0x7fffffff;
	if (ix < 0x3d900000) {	/* |x| < 2**(-38) */
		y = x + c;
		*err = c - (y - x);
		return (y);
	} else {
		z = x * x;
		/* use double precision at a13 and on */
		t = z * a[40];
		for (i = 38; i >= 26; i -= 2)
			t = z * (a[i] + t);
		ee[0] = t;
		x_h = (double) ((float) x);
		z_h = (double) ((float) z);
		x_l = c - (x_h - x);
		z_l = (x_h * x_h - z_h);
		zz[0] = z;
		zz[1] = z_h;
		zz[2] = z_l + x_l * (x + x_h);

		/*
		 * compute (1+z*(a1+z*(a2+z*(a3...+a12+e)))) by call
		 * __libmcr_mx_poly
		 */

		__libmcr_mx_poly(zz, a, ee, 12);

		/* finally x*(1+z*(p1+...)) */
		r = x_h * ee[0];
		t = x * ee[1] + x_l * ee[0];
		s = t + r;
		*err = t - (s - r);
		return (s);
	}
}
