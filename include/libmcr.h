
#pragma ident "@(#)libmcr.h 1.4 04/02/25 SMI"

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
 * This file is for those who want to call the directly the correctly
 * rounded elementary functions as well as the extra precise elementary
 * function with correction term of 20 extra bits.
 */
#ifndef __LIBMCR_H
#define	__LIBMCR_H

#ifdef __cplusplus
extern "C" {
#endif

#pragma ident   "@(#)libmcr.h	1.4    04/02/25"

/*
 * name of the correctly rounded elementary functions
 */
extern double __libmcr_exp(double);
extern double __libmcr_log(double);
extern double __libmcr_pow(double, double);
extern double __libmcr_atan(double);
extern double __libmcr_sin(double);
extern double __libmcr_cos(double);
extern double __libmcr_tan(double);

/*
 * elementary functions with correction term of 20 extra bits
 */
extern double __libmcr_mx_exp(double, double *);
extern double __libmcr_mx_log(double, double *);
extern double __libmcr_mx_pow(double, double, double *);
extern double __libmcr_mx_atan(double, double *);
extern double __libmcr_mx_sin(double, double *);
extern double __libmcr_mx_cos(double, double *);
extern double __libmcr_mx_tan(double, double *);

#ifdef __cplusplus
}
#endif

#endif				/* !defined(__LIBMCR_H) */
