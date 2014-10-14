
#pragma ident "@(#)__libmcr_mm_ln2.c 1.3 04/02/25 SMI"

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
 * __libmcr_mm_ln2 (int *p, int nw)
 *
 * Return p = log2 in nw precision.
 *
 *  Algorithm:
 *	   z+1           1     1      1
 *	ln----- = 2 * [ --- + ---- + ---- + ... ]	        (1)
 *         z-1           z    3z^3   5z^5
 *
 * Set z = 3 in equation (1) we get
 *
 *      ln(2) = 2 * ( 1/3 + 1/(3*3^3) + 1/(5*3^5) + ... )	(2)
 *
 * Since 1/3^2 ~ 0.1111; therefore, every terms in this series will
 * increase about 3.17 bits; kind of slow.
 *
 * A faster way is use (1) to compute  a = log(16/15) (z = 31),
 * b = log(25/24) (z = 49), and  c = log(91/90) (z = 161). We have
 *
 *	a = log(32/30)   =  4*ln2 -   ln3 -   ln5
 *	b = log(50/48)   = -3*ln2 -   ln3 + 2*ln5
 *	c = log(162/160) = -4*ln2 + 4*ln3 -   ln5
 *
 *	or  [ 4  -1  -1 ] [ln2]   [a]
 *	    [-3  -1   2 ]*[ln3] = [b]
 *	    [-4   4  -1 ] [ln5]   [c]
 * Thus
 *	    [ 7   5   3 ] [a]   [ln2]
 *	    [11   8   5 ]*[b] = [ln3]
 *	    [16  12   7 ] [c]   [ln5]
 *
 * and
 *     i    j    k
 * ln(2  * 3  * 5  ) = [i, j, k]* [ ln2 ]
 *				  [ ln3 ]
 *				  [ ln5 ]
 *
 *                   = [i, j, k]*[ 7  5  3] [a]
 *				 [11  8  5]*[b]
 *				 [16 12  7] [c]
 *
 * Thus, as a special case, ln2 = 7a + 5b + 3c.
 *
 * The convergent rate for a,b,c are, respectively, 9.91, 11.23, 14.66 bits.
 *
 * Note: for nw<=40, we use the pre-computed value for ln2, which is up to
 * 912 bits. This is mostly enough for practical use. The faster method
 * is thus not that necessary.
 *
 * (For comparison, to achieve N bit precision,
 * (1) gives N/3.17 which abc gives N(1/9.91+1/11.23+1/14.66). The ratio
 * is 1/3.17 vs 1/9.91+1/11.23+1/14.66  or  0.3155 vs 0.2582; that is, a
 * ~18% speed improvement.)
 *
 * Accuracy:
 *	there is one extra word. The reliable significant bits is 24*(nw-3).
 */

#include <stdlib.h>
#include "mcr.h"
#ifdef DEBUG
#include <stdio.h>
#endif

static const int ln2[] = {	/* for nw <= 40, contain 912 bits of ln2 */
	0x1, -1, 0xB17217, 0xF7D1CF, 0x79ABC9, 0xE3B398,
	0x3F2F6, 0xAF40F3, 0x432672, 0x98B62D, 0x8A0D17, 0x5B8BAA,
	0xFA2BE7, 0xB87620, 0x6DEBAC, 0x985595, 0x52FB4A, 0xFA1B10,
	0xED2EAE, 0x35C138, 0x214427, 0x573B29, 0x1169B8, 0x253E96,
	0xCA1622, 0x4AE8C5, 0x1ACBDA, 0x11317C, 0x387EB9, 0xEA9BC3,
	0xB13660, 0x3B256F, 0xA0EC76, 0x57F74B, 0x72CE87, 0xB19D65,
	0x48CAF5, 0xDFA6BD, 0x383032, 0x48655F,
};

void
__libmcr_mm_ln2(int *p, int nw)
{
	int *ma, *md, *mt, *mp, nwp, i, n;
#ifdef DEBUG
	int j, mtail[4], *mw;
#endif

#ifndef DEBUG
	if (nw <= 40) {
		for (i = 0; i < nw; i++)
			p[i] = ln2[i];
		return;
	}
#endif
	nwp = nw + 1;

	ma = (int *) calloc(nwp, sizeof (int));
	md = (int *) calloc(nwp, sizeof (int));
	mp = (int *) calloc(nwp, sizeof (int));
	mt = (int *) calloc(nwp, sizeof (int));
#ifdef DEBUG
	mw = (int *) calloc(nwp, sizeof (int));
#endif

	/* Set  ma = 1/3, md = 1/9, mp = 1/3 */
	ma[0] = md[0] = mp[0] = 1;
	ma[1] = md[1] = mp[1] = -1;
	for (i = 2; i < nwp; i++) {
		/* ma = 1/3 0.55555555.. in hex */
		mp[i] = ma[i] = 0x555555;
	}
	for (i = 2; i < nwp; i++)
		md[i] = 0x1C71C7;	/* md = 1/9 0.1C71C71C.. in hex */

#ifdef DEBUG
	printf("ma = %1.20e\n", __libmcr_mi_mitod(ma, nw, mtail));
	printf("md = %1.20e\n", __libmcr_mi_mitod(md, nw, mtail));
	printf("mt = %1.20e\n", __libmcr_mi_mitod(mt, nw, mtail));
#endif
	n = 3;
	__libmcr_mm_mul(ma, md, ma, nwp);
	__libmcr_mm_divi(ma, n, mt, nwp);
	__libmcr_mm_add(mt, mp, mp, nwp);
#ifdef DEBUG
	printf("Loop 1:\n");
	printf("ma = %1.20e\n", __libmcr_mi_mitod(ma, nw, mtail));
	printf("md = %1.20e\n", __libmcr_mi_mitod(md, nw, mtail));
	printf("mt = %1.20e\n", __libmcr_mi_mitod(mt, nw, mtail));
	printf("mp = %1.20e\n", __libmcr_mi_mitod(mp, nw, mtail));
#endif
	/* loop until mt is insignificant compare to mp */
	while ((mp[1] - mt[1]) < nwp) {
		n += 2;
		__libmcr_mm_mul(ma, md, ma, nwp);
		__libmcr_mm_divi(ma, n, mt, nwp);
		__libmcr_mm_add(mt, mp, mp, nwp);
#ifdef DEBUG
		__libmcr_mm_add(mp, mp, mw, nwp);
		__libmcr_mm_sub(ln2, mw, mw, 40);
		j = __libmcr_mi_ilogb(ln2) - __libmcr_mi_ilogb(mw);
		printf("At term %2d, %3d bits matched.\n", n / 2, j);
		free(mw);
#endif
	}
	__libmcr_mm_scalbn(mp, nwp, 1);
	for (i = 0; i < nw; i++)
		p[i] = mp[i];

	free(ma);
	free(md);
	free(mp);
	free(mt);
}
