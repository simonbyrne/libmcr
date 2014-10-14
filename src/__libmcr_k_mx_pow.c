
#pragma ident "@(#)__libmcr_k_mx_pow.c 1.7 04/02/25 SMI"

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
 * double __libmcr_k_mx_pow((double*)px,(double*)py,(double*)e, (int*)pexp)
 * Primary routine for the correctly rounded pow function.
 *
 * Method:
 *	1. Compute log(x) = n*ln2 + log(ri) + log(x'/ri)
 *			  = (f1+f2) + (h1+h2) + (g1+g2)
 *
 * 	where n,x',ri are chosen as follow:
 *	    (1) if x>=1, choose n,r such that x = 2**n * r, 1<=r < 2, and
 *		ri = r chopped to 6 significant bits.
 *	    (2) if x<1, choose n,r such that x = 2**n * r, 0.5<=r < 1, and
 *		ri = r rounded up to 6 significant bits.
 * 	(This guarantee the sign of f1,f2, h1,h2, g1,g2 are the same.)
 * 	Implementation details:
 *	if x>=1, n = (hx>>20)-0x3ff, hr = hx-(n<<20), hri = hr&0x7fff8000;
 *	i = (hri-0x3ff00000)>>15, (i=0,...,0x1f=31)
 *	(if i = (hri-0x3fe00000)>>15, then i=32,...,0x1f=63)
 *	if x <1, n=(hx>>20)-0x3fe, hr=hx-(n<<20), hri=(hr+0x7fff)&0x7fff8000;
 *	i = (hri-0x3fe00000)>>15, (i=0,...,0x1f=31)
 *
 *
 *
 *	2. Taylor series for log(r/ri).
 *	   Let s = (r-ri)/(r+ri). Then
 *		log(r/ri) = log((1+s)/(1-s)) = 2s + 2/3 s^3 + 2/5 s^5 + ...
 *	   Since |r-ri|<=0.0315, we have |s|<0.01575, s*s <0.0002480625
 *	   = 2**(-11.977..).
 *	   8 terms is enough for 95 bits accuracy, and one can use double
 *	   precision arithemtic from the fourth term on (the fourth term has
 *	   more than 35 bits shift to the left of s.)
 *
 *		     3        5        7        9         11        13        15
 *	   2s + 2/3 s  + 2/5 s  + 2/7 s  + 2/9 s  + 2/11 s  + 2/13 s  + 2/15 s
 *			          ^^^^^^^^^^^^^^^^
 *			      from this terms on use double precision arithmetic
 *	   Note:
 *		How to compute s1 + s2 = s, where s1 is the first 35 bits of s.
 *		Let p = r-ri (exact), ra = first 17 bits of r, rb = second 18
 *		bits of r, rc = last 18 bits of r; then q1 = ra+ri is exact and
 *		has at most 18 bits, and (r+ri)=q1+rb+rc.
 *		Let qi = 1/(r+ri). Then
 *			s = p*qi; s1 = s rounded to 35 bit;
 *			s2 = p/(q1+rb+rc)-s1 = qi*(p-s1*q1-s1*rb-s1*rc)
 *			   = qi*(((p-s1*q1)-s1*rb)-s1*rc)
 *		How to compute w1+w2= (s1+s2)**2
 *			w1 = s*s rounded to 35 bits
 *			s11 = first 17 bits of s1
 *			s12 = last 18 bits of s1
 *			w2 = (s1*(s11+s12)+s1*s2+s2*s)-w1
 *			   = (((s1*s11-w1)+s1*s12)+s1*s2)+s2*s:
 *
 *		Then log(r/ri) is computed by (let w = s*s rounded)
 *		    1. t = w*(2/7+w*(...2/15))	... in double precision
 *		    2. compute (q1+q2) = w*(2/5+t):
 *			p  = 2/5+t;
 *			p1 = p chopped to 35 bits
 *			p2 = (2/5-p1)+t
 *			q  = w*p
 *			q1 = q chopped to 35 bits
 *			q2 = (w1+w2)*(p1+p2)-q1
 *			   = ((w11*p1-q1)+w12*p1)+w1*p2+w2*p
 *		    3. compute (q1+q2) = w*(2/3+(q1+q2))
 *			p  = 2/3+q;
 *			p1 = p chopped to 35 bits
 *			p2 = ((2/3-p1)+q1)+q2;
 *			q  = w*p
 *			q1 = q chopped to 35 bits
 *			q2 = (w1+w2)*(p1+p2)-q1
 *			   = ((w11*p1-q1)+w12*p1)+w1*p2+w2*p
 *		    4. compute (g1+g2) = s*(2.0+(q1+q2))
 *			p  = 2.0+q;
 *			p1 = p chopped to 35 bits
 *			p2 = ((2-p1)+q1)+q2;
 *			q  = s*p
 *			g1 = q
 *			g2 = (s1+s2)*(p1+p2)-g1
 *			   = ((s11*p1-g1)+s12*p1)+s1*p2+s2*p
 *
 *	3. Table look up for log(ri)
 *		if x>=1, n = (hx>>20)-0x3ff, hr = hx-(n<<20),
 *		hri = hr&0x7fff8000; i = (hri-0x3ff00000)>>15, (i=0,...,0x1f=31)
 *		if x <1, n = (hx>>20)-0x3fe, hr = hx-(n<<20),
 *		hri = (hr+0x7fff)&0x7fff8000;
 *		i = 0x20|((hri-0x3fe00000)>>15), (i=32,...,0x3f=63)
 *
 *		h1 = lnri[2i]
 *		h2 = lnri[2i+1]
 *
 *	4. n*ln(2)
 *	   since n at most 11 bits, let ln2_h = ln2 chopped to 42 bits,
 *	   ln2_l = (ln2-ln2_h). Then we have ~42+53 = 95 accuracy.
 *		n*ln2 = f1 + f2, where f1=n*ln2_h, f2 = n*ln2_l
 *
 *
 *	5. Sum them together
 *		[      f1      ] [      f2      ]
 *		     [     h1    ] [      h2     ]
 *			 [   g1     ] [     g2     ]
 *		--------------------------------------------------
 *			z1             z2
 *
 * 	Because of all f1,f2,g1,g2,h1,h2 are the same sign, we can add them by
 *	    (1) if n!=0,
 *		t = f2+(h2+g2), z1 = f1+(h1+(g1+t)); z2 = (((f1-z1)+h1)+g1)+t);
 *	    (2) if n=0,(i&0x1f)!=0,
 *		t = (h2+g2), z1 = (h1+g1)+t; z2 = ((h1-z1)+g1)+t;
 *	    (3) otherwise
 *		z1 = g1+g2; z2=(g1-z1)+g2;
 *
 */


static const double lnri[] = {	/* (g1,g2) = log(ri) */
	-6.93147180559945286227e-01, -2.31904681384629955842e-17,
	-6.62375521893191598899e-01, -2.21472949355623993815e-17,
	-6.32522558743510376900e-01, -8.99370045836625817055e-17,
	-6.03535021870258092669e-01, -8.41284323033855253542e-17,
	-5.75364144903561802735e-01, -5.21432123288512733580e-17,
	-5.47965170715447369432e-01, -4.27036249710694290829e-17,
	-5.21296923633286057864e-01, -2.92129219594743645195e-17,
	-4.95321437230025418685e-01, -1.03692737654828536582e-17,
	-4.70003629245735521369e-01, -3.22817387357877884699e-17,
	-4.45311016655364044770e-01, -7.86710210153660525541e-18,
	-4.21213465076303528178e-01, -2.24071485007655527051e-17,
	-3.97682967666109388194e-01, -4.48365767425228952480e-17,
	-3.74693449441410642020e-01, -5.15868400023945866421e-17,
	-3.52220593589352093389e-01, -5.72333169491824846705e-18,
	-3.30241686870576811597e-01, -4.46828295937739644955e-17,
	-3.08735481649613230370e-01, -3.93119651461097186189e-17,
	-2.87682072451780901368e-01, -2.60716061644256366790e-17,
	-2.67062785249045198110e-01, -4.81822359039376529103e-17,
	-2.46860077931525784267e-01, -1.36174337174836786303e-17,
	-2.27057450635346075307e-01, -9.55141576273848843149e-18,
	-2.07639364778244489562e-01, -1.20532432166861289490e-17,
	-1.88591169807550002036e-01, -2.03234113964319869169e-17,
	-1.69899036795397445632e-01, -2.72687747391850040238e-17,
	-1.51549898127200932674e-01, -5.16695936846155866966e-18,
	-1.33531392624522599055e-01, -2.40911179519688264552e-17,
	-1.15831815525121700761e-01, -4.33848436980809518520e-18,
	-9.84400728132525104641e-02, -9.43877817413932010836e-18,
	-8.13456394539524008103e-02, -5.07707635593116915980e-18,
	-6.45385211375711642656e-02, -7.40730114612152221714e-18,
	-4.80092191863606063129e-02, -1.43909033472922027777e-18,
	-3.17486983145802981188e-02, -3.03822630846808578765e-18,
	-1.57483569681391676054e-02, -1.00215786305289717445e-18,
	0.00000000000000000000e+00, 0.00000000000000000000e+00,
	3.07716586667536873279e-02, 1.04317320290059658791e-18,
	6.06246218164348399382e-02, 2.64240259387269303297e-18,
	8.96121586896871241690e-02, 8.45097487414974168052e-18,
	1.17783035656383441858e-01, 1.26806192330550878698e-17,
	1.45182009844497889040e-01, 8.24241878302247384888e-18,
	1.71850256926659200607e-01, 2.17331217946175414938e-17,
	1.97825743329919867541e-01, 1.28211943729801403852e-17,
	2.23143551314209737102e-01, 1.86643050183041144619e-17,
	2.47836163904581241457e-01, 1.53233660369263887880e-17,
	2.71933715483641758048e-01, 7.83319637697442012432e-19,
	2.95464212893835842522e-01, 3.38650426271979242758e-17,
	3.18453731118534588695e-01, 2.71147793673262328817e-17,
	3.40926586970593192838e-01, 1.74671364435447471171e-17,
	3.62905493689368419119e-01, 3.40187897759468488654e-17,
	3.84411698910332000345e-01, 3.93896542236110947420e-17,
	4.05465108108164329348e-01, 5.26300132052951828449e-17,
	4.26084395310900032605e-01, 3.05193834657831604506e-17,
	4.46287102628419474204e-01, 3.73286100366082289238e-17,
	4.66089729924599183164e-01, 4.13946279913534206633e-17,
	4.85507815781700768909e-01, 3.88928005374057755235e-17,
	5.04556010752395200925e-01, 8.61337835889177368731e-17,
	5.23248143764547757328e-01, 7.91884202461647228477e-17,
	5.41597282432744298042e-01, 7.35346600012592639357e-17,
	5.59615787935422659416e-01, 2.68549258021230795581e-17,
	5.77315365034823502199e-01, 1.02118710615541637849e-16,
	5.94707107746692775763e-01, 1.37516899643236739351e-17,
	6.11801541105992829905e-01, 7.36245430137896526752e-17,
	6.28608659422374094206e-01, 4.35387426079703807147e-17,
	6.45137961373584700731e-01, 9.34696092012090558654e-19,
	6.61398482245364904841e-01, 1.03418968676881641854e-16,
	6.77398823591806031885e-01, 1.08924484074250449664e-16,
};

static const double
ln2_h = 6.93147180559890330187e-01,	/* 3FE62E42 FEFA3800 */
ln2_l = 5.49792301870837115524e-14,	/* 3D2EF357 93C76730 */
v1 = 6.66666666666666629659e-01,	/* 2/3  3FE55555 55555555 */
v1t = 3.70074341541718826264e-17,	/* 3C855555 55555555 */
v2 = 3.99999999999999966693e-01,	/* 2/5  3FD99999 99999999 */
v2t = 3.33066907387546949801e-17,	/* 3C833333 33333333 */
v3 = 2.85714285714285698425e-01,	/* 2/7  FD24924 92492492 */
v4 = 2.22222222222222209886e-01,	/* 2/9  3FCC71C7 1C71C71C */
v5 = 1.81818181818181795472e-01,	/* 2/11 3FC745D1 745D1745 */
v6 = 1.53846153846153826938e-01,	/* 2/13 3FC3B13B 13B13B13 */
v7 = 1.33333333333333331482e-01,	/* 2/15 3FC11111 11111111 */
two54 = 18014398509481984.0, one = 1.0, zero = 0.0;

#include "mcr.h"

#ifdef DEBUG
#include <stdio.h>
#define	H4(x)  *(int *)&x, *(1+(int *)&x), *(2+(int *)&x), *(3+(int *)&x)
#define	H2(x)  *(int *)&x, *(1+(int *)&x)
#endif

double
__libmcr_k_mx_pow(double *px, double *py, double *err, int *pexp)
{
	int hr, hri, j, n;
	double x = *px, y = *py;
	double r, h1, h2, ri, y1, y2, p, q, ra, rb, rc, q1, qi, q2, s,
	    v, w12, w11;
	double f1, f2, u, w, s1, s2, w1, w2, z1, z2, t1, t2, g1, g2, t,
	    s11, s12, p1, p2;
#ifdef DEBUG
	long double xx, yy, zz, ww, logl(), expl(), scalbnl();
	double log2();
#endif

	r = x;
	hr = HIGH_WORD(r);
#if DEBUG
	printf("\nx = %08X %08X, log2(x) = %g\n", H2(x), log2(x));
	printf("y = %08X %08X, %1.20e\n", H2(y), y);
#endif
	*err = 0.0;
	*pexp = 0;
	/* filter out exception argument */

	/* compute n,x,ri */
	if ((n = (hr >> 20) - 0x3ff) >= 0) {
		hr = hr - (n << 20), hri = hr & 0x7fff8000;
	} else {
		if ((n + 0x3ff) == 0) {	/* subnormal */
			r *= two54;
			hr = HIGH_WORD(r);
			n = (hr >> 20) - 0x3fe;
			hr = hr - (n << 20);
			n -= 54;
		} else {
			n += 1;
			hr = hr - (n << 20);
		}
		hri = (hr + 0x7fff) & 0x7fff8000;	/* rounded up */
	}
	ri = zero;
	HIGH_WORD(r) = hr;
	HIGH_WORD(ri) = hri;
	qi = one / (r + ri);
	p = r - ri;
	ra = zero;
	HIGH_WORD(ra) = hr & 0xfffffff0;	/* ra = first 17 bits of r */
	rb = r;
	q1 = ri + ra;
	LOW_WORD(rb) = LOW_WORD(x) & 0xfffc0000; /* last 18 bits of r is zero */
	rc = r - rb;
	rb = rb - ra;

	t = (double) n;
	j = (hri - 0x3fe00000) >> 14;
	h1 = lnri[j];
	h2 = lnri[j + 1];
	f1 = t * ln2_h;
	f2 = t * ln2_l;

#if DEBUG
	printf("n    = %d\n", n);
	printf("j    = %d\n", j);
	printf("r    = %08X %08X\n", H2(r));
	printf("q1   = %08X %08X\n", H2(q1));
	printf("ri   = %08X %08X\n", H2(ri));
	printf("ra   = %08X %08X\n", H2(ra));
	printf("rb   = %08X %08X\n", H2(rb));
	printf("rc   = %08X %08X\n", H2(rc));
	v = ra + rb + rc;
	printf("rabc = %08X %08X\n", H2(v));
	xx = ri;
	zz = logl(xx);
	ww = (long double) h1 + (long double) h2;
	printf("h1,h2    = %08X %08X %08X %08X\n", H2(h1), H2(h2));
	printf("h1+h2    = %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));

	xx = n;
	zz = xx * logl(2.0L);
	ww = (long double) f1 + (long double) f2;
	printf("f1,f2    = %08X %08X %08X %08X\n", H2(f1), H2(f2));
	printf("f1+f2    = %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif
	/* compute s = (p/(q1+rb+rc) - s1 */
	s11 = s1 = s = p * qi;
	LOW_WORD(s1) &= 0xfffc0000;	/* s1 = 35 bits of s */
	t1 = s1 * q1;
	t2 = s1 * rb;
	u = (p - t1) - t2;
	s2 = qi * (u - s1 * rc);
#if DEBUG
	printf("p  = %08X %08X %1.20e\n", p, p);
	printf("qi = %08X %08X %1.20e\n", qi, qi);
	printf("q1 = %08X %08X %1.20e\n", q1, q1);
	printf("s  = %08X %08X %1.20e\n", s, s);
	printf("s1 = %08X %08X %1.20e\n", s1, s1);
	printf("s2 = %08X %08X %1.20e\n", s2, s2);
	xx = r;
	yy = ri;
	zz = (xx - yy) / (xx + yy);
	ww = (long double) s1 + (long double) s2;
	printf("s1,s2    = %08X %08X %08X %08X\n", H2(s1), H2(s2));
	printf("s1+s2    = %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif

	/* compute w = s*s */
	w1 = w = s * s;
	LOW_WORD(w1) &= 0xfffc0000;	/* w1 = 35 bits of w */
	HIGH_WORD(s11) &= 0xfffffff0;
	LOW_WORD(s11) = 0;	/* s11 = first 17 bits of s1 */
	s12 = s1 - s11;
	t1 = s1 * s11;
	u = s1 * s12 - (w1 - t1);
	v = s1 * s2 + s2 * s;
	w2 = u + v;
#if DEBUG
	printf("s11= %08X %08X %1.20e\n", s11, s11);
	printf("s12= %08X %08X %1.20e\n", s12, s12);
	t = s11 + s12;
	printf("s12+s11= %08X %08X %1.20e\n", t, t);
	printf("s1     = %08X %08X %1.20e\n", s1, s1);
	zz = zz * zz;
	ww = (long double) w1 + (long double) w2;
	printf("w1,w2    = %08X %08X %08X %08X\n", H2(w1), H2(w2));
	printf("w1+w2    = %08X %08X %08X %08X\n", H4(ww));
	ww = (long double) s *(long double) s;
	printf("s*s      = %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif

	/* compute t = w*(2/7+w*(...2/15))	... in double precision */
	t = w * (v3 + w * (v4 + w * (v5 + w * (v6 + w * v7))));
	w11 = zero;
	HIGH_WORD(w11) = HIGH_WORD(w1) & 0xfffffff0;
	w12 = w1 - w11;

	/* compute (q1+q2) = w*(2/5+t) */
#if DEBUG
	yy = 2.0L / 5.0L;
	zz = (long double) v2 + (long double) v2t;
	printf("  (2/5  )= %08X %08X %08X %08X\n", H4(zz));
	printf("expected = %08X %08X %08X %08X\n", H4(yy));
#endif
	p1 = p = v2 + t;
	q1 = q = w * p;
	LOW_WORD(p1) &= 0xfffc0000;
	p2 = (t - (p1 - v2)) + v2t;
#if DEBUG
	ww = (long double) w1 + (long double) w2;
	yy = (long double) t + 2.0L / 5.0L;
	zz = (long double) p1 + (long double) p2;
	printf("  (2/5+t)= %08X %08X %08X %08X\n", H4(zz));
	printf("expected = %08X %08X %08X %08X\n", H4(yy));
#endif
	LOW_WORD(q1) &= 0xfffc0000;
	u = w2 * p + w1 * p2;
	v = w12 * p1 - (q1 - w11 * p1);
	q2 = u + v;
#if DEBUG
	yy = (long double) t + 2.0L / 5.0L;
	yy = ww * yy;
	zz = (long double) q1 + (long double) q2;
	printf("w*(2/5+t)= %08X %08X %08X %08X\n", H4(zz));
	printf("expected = %08X %08X %08X %08X\n", H4(yy));
#endif

	/* compute (q1+q2) = w*(2/3+(q1+q2)) */
	q2 += v1t;
	p1 = p = v1 + q;
	LOW_WORD(p1) &= 0xfffc0000;
	p2 = q1 - (p1 - v1);
	p2 += q2;
	q1 = q = w * p;
	LOW_WORD(q1) &= 0xfffc0000;
	u = w2 * p + w1 * p2;
	v = w12 * p1 - (q1 - w11 * p1);
	q2 = u + v;

#if DEBUG
	yy += 2.0L / 3.0L;
	yy = ww * yy;
	zz = (long double) q1 + (long double) q2;
	printf("w*(2/3+.)= %08X %08X %08X %08X\n", H4(zz));
	printf("expected = %08X %08X %08X %08X\n", H4(yy));
#endif
	/* compute (q1+q2) = s*(2.0+(q1+q2)) */
	p1 = p = 2.0 + q;
	LOW_WORD(p1) &= 0xfffc0000;
	p2 = q1 - (p1 - 2.0);
	p2 += q2;
#if DEBUG
	yy += 2.0L;
	zz = (long double) p1 + (long double) p2;
	printf("(2 + .)  = %08X %08X %08X %08X\n", H4(zz));
	printf("expected = %08X %08X %08X %08X\n", H4(yy));
#endif
	q = s * p;
	g1 = q;
	u = s1 * p2 + s2 * p;
	g2 = (s12 * p1 - (g1 - s11 * p1)) + u;

#if DEBUG
	yy = ((long double) s1 + (long double) s2) * yy;
	zz = (long double) g1 + (long double) g2;
	printf("s*(2 + .)= %08X %08X %08X %08X\n", H4(zz));
	printf("expected = %08X %08X %08X %08X\n", H4(yy));
#endif

#if DEBUG
	xx = r;
	yy = ri;
	zz = logl(xx / yy);
	ww = (long double) g1 + (long double) g2;
	printf("log(r/ri)= %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));
	ww = (long double) (s *
		(2.0 + w *
		(v1 + w *
		(v2 + w * (v3 + w * (v4 + w * (v5 + w * (v6 + w * v7))))))));
	printf("formula  = %08X %08X %08X %08X\n", H4(ww));
#endif
	/* now add (z1,z2) = (f1+f2)+(h1+h2)+(g1+g2) */

#if DEBUG
	zz = (long double) f1 + (long double) f2 + (long double) h1 +
		(long double) h2 + (long double) g1 + (long double) g2;
#endif
	if (n == 0) {
		if ((j & 0x3f) == 0) {
			z1 = g1 + g2;
			LOW_WORD(z1) &= 0xfffc0000;
			z2 = g2 - (z1 - g1);
#if DEBUG
			ww = (long double) z1 + (long double) z2;
			printf("    g    = %08X %08X %08X %08X\n", H4(ww));
			printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif
		} else {
			t = (h2 + g2);
			z1 = h1 + g1 + t;
			LOW_WORD(z1) &= 0xfffc0000;
			z2 = g1 - (z1 - h1);
			z2 += t;
#if DEBUG
			ww = (long double) z1 + (long double) z2;
			printf("h + g    = %08X %08X %08X %08X\n", H4(ww));
			printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif
		}
	} else {
		t = f2 + (h2 + g2);
		z1 = f1 + h1 + g1 + t;
		LOW_WORD(z1) &= 0xfffc0000;
		z2 = h1 - (z1 - f1);
		z2 += g1;
		z2 += t;
#if DEBUG
		ww = (long double) z1 + (long double) z2;
		printf("f+h+g    = %08X %08X %08X %08X\n", H4(ww));
		printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif
	}
#if DEBUG
	zz = logl((long double) x);
	ww = (long double) z1 + (long double) z2;
	printf("log(x)   = %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif

	/* compute y*(z1+z2) */
	t = y * z2;
	y1 = zero;
	HIGH_WORD(y1) = HIGH_WORD(y) & 0xfffffff0;
	y2 = y;
	LOW_WORD(y2) &= 0xfffc0000;
	t += (y - y2) * z1;
	t1 = z1 * y1;
	t2 = z1 * (y2 - y1);
	w1 = (t1 + (t2 + t));
	w2 = (t2 - (w1 - t1)) + t;
#if DEBUG
	printf("y = %08X %08X %1.20e\n", H2(y), y);
	printf("z1= %08X %08X %1.20e\n", H2(z1), z1);
	printf("z2= %08X %08X %1.20e\n", H2(z2), z2);
	printf("t = %08X %08X %1.20e\n", H2(t), t);
	printf("t1= %08X %08X %1.20e\n", H2(t1), t1);
	printf("t2= %08X %08X %1.20e\n", H2(t2), t2);
	printf("w1= %08X %08X %1.20e\n", H2(w1), w1);
	printf("w2= %08X %08X %1.20e\n", H2(w2), w2);
	zz = (long double) y *zz;
	ww = (long double) w1 + (long double) w2;
	printf("y*log(x) = %08X %08X %08X %08X\n", H4(ww));
	printf("expected = %08X %08X %08X %08X\n", H4(zz));
#endif

	/* call __libmcr_k_mx_exp */
	t1 = __libmcr_k_mx_exp(w1, &q, pexp);
#if DEBUG
	printf("(t1,t2),exponent: %08X %08X %08X %08X %d\n",
		H2(t1), H2(q), *pexp);
#endif

	/* adjustment: (t1+q)*exp(w2)=(t1+q)(1+w2)=t1+q+t1*w2 */
	if (q == zero) {
		p = t1;
		*err = zero;
	} else {
		q += t1 * w2;
		p = t1 + q;
		*err = q - (p - t1);
	}
#if DEBUG
	q = *err;
	printf("p = __libmcr_k_mx_exp  %08X %08X %1.20e\n", H2(p), p);
	printf("q = err                %08X %08X %1.20e\n", H2(q), q);
	printf("exponent    %d\n", *pexp);
	zz = scalbnl(expl((long double) ww), -(*pexp));
	ww = (long double) p + (long double) q;
	printf("p, q     = %08X %08X %08X %08X\n", H2(p), H2(q));
	printf("p +q     = %08X %08X %08X %08X %1.30Le\n", H4(ww), ww);
	printf("expected = %08X %08X %08X %08X %1.30Le\n", H4(zz), zz);
#endif

	return (p);

}
