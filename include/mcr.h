
#pragma ident "@(#)mcr.h 1.8 04/02/25 SMI"

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

#pragma ident "@(#)mcr.h	1.8    04/02/25"
/* INDENT OFF */

#include <math.h>	/* fabs, sqrt, sin, cos, tan, ... */
#include "libmcr.h"

/* internal multi-precision routines */
extern void __libmcr_mm_sub(int *, int *, int *, int);
extern void __libmcr_mm_add(int *, int *, int *, int);
extern void __libmcr_mm_div(int *, int *, int *, int);
extern void __libmcr_mm_divi(int *, int,  int *, int);
extern void __libmcr_mm_mul(int *, int *, int *, int);
extern void __libmcr_mm_muli(int *, int,  int *, int);
extern void __libmcr_mm_scalbn(int *, int, int);
extern void __libmcr_mm_sqrt(int *, int *, int);
extern void __libmcr_mm_atan(int *, int *, int);
extern void __libmcr_mm_pio2(int *, int);
extern void __libmcr_mm_pio4(int *, int);
extern void __libmcr_mm_ln2(int *, int);
extern void __libmcr_k_mm_exp(int *, int *, int);
extern void __libmcr_k_mm_sin(int *, int *, int);
extern void __libmcr_k_mm_cos(int *, int *, int);
extern void __libmcr_k_mi_log(double, int *, int);
extern double __libmcr_mi_atan(double, int, int *, int);
extern double __libmcr_mi_exp(double, int, int *, int);
extern double __libmcr_mi_pow(double, double, int, int *, int);
extern double __libmcr_mi_log(double, int, int *, int);
extern double __libmcr_mi_sin(double, int, int *, int);
extern double __libmcr_mi_cos(double, int, int *, int);
extern double __libmcr_mi_tan(double, int, int *, int);
extern int  __libmcr_k_mi_rem_pio2(double, int *, int);

/* internal 20 bits extra-precise double routines */
extern double __libmcr_k_mx_exp(double, double *, int *);
extern double __libmcr_k_mx_pow(double *, double *, double *, int *);
extern double __libmcr_k_mx_sin(double, double, double *);
extern double __libmcr_k_mx_cos(double, double, double *);
extern double __libmcr_k_mx_tan(double, double, double *);

/* internal auxilary routines */
extern int    __libmcr_mx_check(double *, int, double);
extern int    __libmcr_k_mx_rem_pio2(double, double *);
extern int    __libmcr_k_rem_pio2m(double *, double *, int, int, int);
extern int    __libmcr_mi_ilogb(int *);
extern void   __libmcr_mx_poly(double *, double *, double *, int);
extern void   __libmcr_mi_dtomi(double, int *, int);
extern void   __libmcr_mi_itomi(int, int *, int);
extern void   __libmcr_mi_format(int *, int);
extern double __libmcr_mi_mitod(int *, int, int *);
extern double __libmcr_mi_final(int *, int, int *, int);
extern double __libmcr_k_exactsqrt(double, int *);

/* internal data constant */
extern const double __libmcr_TBL_atan_hi[];
extern const double __libmcr_TBL_atan_lo[];
extern const int __libmcr_TBL_pio2[];
extern const int __libmcr_TBL_pio4[];
extern const int __libmcr_TBL_ipio2[];

/* fsr definitions */
#define	RDmask	0xc0000000
#define	RDrm	0xc0000000
#define	RDrp	0x80000000
#define	RDrz	0x40000000
#define	RDrn	0x00000000

#define	TEmask	0x0f800000
#define	TEnvm	0x08000000
#define	TEofm	0x04000000
#define	TEufm	0x02000000
#define	TEdzm	0x01000000
#define	TEnxm	0x00800000

#define	NS	0x00400000

#define	EXmask	0x000003e0
#define	EXnva	0x00000200
#define	EXofa	0x00000100
#define	EXufa	0x00000080
#define	EXdza	0x00000040
#define	EXnxa	0x00000020

/* macros to extract HIGH and LOW words from a double */
#if defined(__LITTLE_ENDIAN) || defined(__i386)
#define	HIGH_WORD(x) (*(1+(int *)&(x)))
#define	LOW_WORD(x) (*(unsigned *)&(x))
#else
#define	HIGH_WORD(x) (*(int *)&(x))
#define	LOW_WORD(x) (*(1+(unsigned *)&(x)))
#endif
