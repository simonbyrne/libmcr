
#pragma ident "@(#)__libmcr_mi_mitod.c 1.4 04/02/25 SMI"

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
 * double __libmcr_mi_mitod((int*)mx,(int)nw,(int*)mtail)
 * Truncate a mi number to double with guard,sticky,ms0,ms1 in
 * mtail[] to help determine the rounded result.
 *
 * According to our analysis, all basic routines in mm*.c produce results
 * with 24*(nw-3) significant bits. It is thus safe to allow UNC(=10) bits of
 * uncertainty for routine in __libmcr_mi_*.c.
 *
 * Four values will be given by the subroutine __libmcr_mi_mitod in array
 * mtail[]:
 *
 *    |---------------------------------|
 *    | 1.xxx ....xx   r   ss........s  |
 *    |---------------------------------|
 *
 *     bit        bit bit bit       bit
 *      0         53  54  55 -> 24(nw-3)-UNC
 *                   (mrd) (mst,ms0,ms1)
 * mtail[0] = mrd --- round bit (54th bit)
 * mtail[1] = mst --- sticky bit (union of the 55th bit to the last bits)
 * mtail[2] = ms0 --- =1 if all s are 0
 * mtail[3] = ms1 --- =1 if all s are 1
 *
 * Note: if 24(nw-3)-UNC < 55, then we set ms0=ms1=1.
 *
 * We assume UNC < 19 in the following algorithm.
 */

#include "mcr.h"
#ifdef DEBUG
#include <stdio.h>
#endif

#define	UNC  10

static const int m24 = 0xffffff;

double
__libmcr_mi_mitod(mx, nw, mtail)
	int mx[], nw, mtail[];
{
	double x;
	int iexp, i, j, k, hx, lx, m2, m3, m4, m5, m6, m7, mgd, mst, ms0,
	    ms1, m2b;

	/* clean up mx */
	__libmcr_mi_format(mx, nw);

	m2 = mx[2];
#ifdef DEBUG
	printf("\n");
#endif
	/* check for zero,over/underflow cases */
	if (m2 == 0 || m2 >= 0x7ff00000) {	/* 0,inf or nan */
		HIGH_WORD(x) = (mx[0] & 0x80000000) | m2;
		LOW_WORD(x) = 0;
		mtail[0] = mtail[1] = mtail[2] = mtail[3] = 0;
		return (x);
	}
	i = mx[1];
	if (i < -45) {		/* underflow */
		HIGH_WORD(x) = (mx[0] & 0x80000000);
		LOW_WORD(x) = 0;
		mtail[0] = mtail[2] = mtail[3] = 0;
		mtail[1] = 1;
		return (x);
	}
	if (i > 42) {		/* overflow */
		HIGH_WORD(x) = (mx[0] & 0x80000000) | 0x7ff00000;
		LOW_WORD(x) = 0;
		mtail[0] = mtail[1] = mtail[2] = mtail[3] = 0;
		return (x);
	}
	/* compute m2b = # of significant bits in m2 */
	m2b = 0;
	j = 1;
	while (j <= m2) {
		m2b += 1;
		j <<= 1;
	}

#ifdef DEBUG
	printf("m2b = %d\n", m2b);
#endif
	ms0 = 1;
	ms1 = 1;

	/* compute ms0,ms1 from mx[6] to mx[nw-3] when nw > 8 */
	for (i = 6; i < (nw - 2) && (ms0 + ms1) != 0; i++) {
		if (mx[i] != m24)
			ms1 = 0;
		if (mx[i] != 0)
			ms0 = 0;
	}
	mst = (1 - ms0);

	/*
	 * put the remaining mx to m2-m7, where m6 and m7 contain the last
	 * two term of mx
	 */
	if (nw >= 8) {
		m3 = mx[3];
		m4 = mx[4];
		m5 = mx[5];
		m6 = mx[nw - 2];
		m7 = mx[nw - 1];
	} else {
		m7 = 0;
		switch (nw) {
		case 7:
			m3 = mx[3];
			m4 = mx[4];
			m5 = mx[5];
			m6 = mx[6];
			break;
		case 6:
			m3 = mx[3];
			m4 = mx[4];
			m5 = mx[5];
			m6 = 0;
			break;
		case 5:
			m3 = mx[3];
			m4 = mx[4];
			m5 = m6 = 0;
			break;
		case 4:
			m3 = mx[3];
			m4 = m5 = m6 = 0;
			break;
		default:
		case 3:
			m3 = m4 = m5 = m6 = 0;
			break;
		}
	}
	mst |= (m6 | m7);

#ifdef DEBUG
	printf("m2,m3,...,  m7 = %06X %06X %06X %06X %06X %06X\n",
		m2, m3, m4, m5, m6, m7);
#endif

	/* normalize to shift m2 to 1 bit */
	k = m2b - 1;
	if (k > 0) {
		i = 32 - k;
		m7 = ((unsigned) (m6 << i) >> 8) | (m7 >> k);
		m6 = ((unsigned) (m5 << i) >> 8) | (m6 >> k);
		m5 = ((unsigned) (m4 << i) >> 8) | (m5 >> k);
		m4 = ((unsigned) (m3 << i) >> 8) | (m4 >> k);
		m3 = ((unsigned) (m2 << i) >> 8) | (m3 >> k);
		m2 >>= k;
	};

#ifdef DEBUG
	printf("after shifting : %06X %06X %06X %06X %06X %06X\n",
		m2, m3, m4, m5, m6, m7);
#endif
	/* exponent */
	i = mx[1];
	iexp = (i + i + i) << 3; /* i*24 */
	iexp += k;
	if (iexp >= 1024) {	/* overflow */
		HIGH_WORD(x) = (mx[0] & 0x80000000) | 0x7ff00000;
		LOW_WORD(x) = 0;
		mtail[0] = mtail[1] = mtail[2] = mtail[3] = 0;
		return (x);
	}
#ifdef DEBUG
	printf("mx[1] = %d\n", mx[1]);
	printf("iexp  = %d\n", iexp);
#endif

	/* guard bit and sticky word */
	j = 0x7ffff;
	mgd = (m5 >> 19) & 1;
	mst |= (m5 & j);

	/* mantissa */
	hx = m3 >> 4;
	lx = (m3 << 28) | (m4 << 4) | ((unsigned) m5 >> 20);

#ifdef DEBUG
	printf("Mantissa: %08X %08X\n", hx, lx);
#endif

#ifdef DEBUG
	printf("ms0,ms1 = %d %d\n", ms0, ms1);
#endif
	/* ms0 and ms1 */
	m5 &= j;
	if ((ms0 + ms1) != 0) {
		i = m24 >> UNC;	/* remember, UNC < 19 */
		if (nw >= 8) {
			k = m7 >> UNC;
			if (m5 != j || m6 != m24 || k != i)
				ms1 = 0;
			if ((m5 | m6 | k) != 0)
				ms0 = 0;
		} else if (nw == 7) {
			k = m6 >> UNC;
			if (m5 != j || k != i)
				ms1 = 0;
			if ((m5 | k) != 0)
				ms0 = 0;
		} else if (nw == 6) {
			j >>= UNC;
			k = m5 >> UNC;
			if (k != j)
				ms1 = 0;
			if (k != 0)
				ms0 = 0;
#ifdef DEBUG
			printf("nw=6: m5=%06X  i=m24>>UNC, k=m5>>UNC=%06X\n",
				m5, k);
			printf("ms0,ms1 = %d %d\n", ms0, ms1);
#endif
		}
		/* for nw <= 5, ms0 and ms1 will remain equal to 1 */
	}
	/* subnormal or underflow (chopped result) */
	if (iexp < -1022) {
#ifdef DEBUG
		printf("subnormal/underflow:\n");
		printf("iexp  = %d\n", iexp);
#endif

		iexp = -1022 - iexp;
		if (iexp >= 53) {	/* underflow */
			HIGH_WORD(x) = mx[0] & 0x80000000;
			LOW_WORD(x) = 0;
			if (iexp == 53) {
				mgd = 1;
				mst |= hx | lx | mgd;
			} else {
				mgd = 0;
				mst = 1;
			}
		} else {
			hx |= 0x100000;
			if (iexp > 20) {
				HIGH_WORD(x) = mx[0] & 0x80000000;
				if (iexp > 32) {
					if (lx != 0)
						ms0 = 0;
					if (lx != -1)
						ms1 = 0;
					j = hx << (65 - iexp);
					mst = j | lx | mgd | mst;
					mgd = (hx >> (iexp - 33)) & 1;
					lx = hx >> (iexp - 32);
					if (j != 0)
						ms0 = 0;
					if ((j >> (65 - iexp)) != -1)
						ms1 = 0;
				} else {
					j = lx << (33 - iexp);
					mst = j | mgd | mst;
					mgd = (lx >> (iexp - 1)) & 1;
					if (iexp == 32)
						lx = hx;
					else
						lx = ((unsigned) lx >> iexp) |
							(hx << (32 - iexp));
					if (j != 0)
						ms0 = 0;
					if ((j >> (33 - iexp)) != -1)
						ms1 = 0;
				}
				LOW_WORD(x) = lx;
			} else {
				j = lx << (32 - iexp);
				mst |= mgd;
				mgd = ((unsigned) j) >> 31;
				j <<= 1;
				mst |= j;
				if (j != 0)
					ms0 = 0;
				if ((j >> (33 - iexp)) != -1)
					ms1 = 0;
				lx = ((unsigned) lx >> iexp) |
					(hx << (32 - iexp));
				hx >>= iexp;
				LOW_WORD(x) = lx;
				HIGH_WORD(x) = (mx[0] & 0x80000000) | hx;
			}
		}
		mtail[0] = mgd;
		mtail[1] = mst;
		mtail[2] = ms0;
		mtail[3] = ms1;
		return (x);
	}
	hx = (mx[0] & 0x80000000) | (hx + ((iexp + 0x3ff) << 20));
	LOW_WORD(x) = lx;
	HIGH_WORD(x) = hx;
	mtail[0] = mgd;
	mtail[1] = mst;
	mtail[2] = ms0;
	mtail[3] = ms1;
	return (x);
}
