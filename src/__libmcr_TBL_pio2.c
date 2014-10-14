
#pragma ident "@(#)__libmcr_TBL_pio2.c 1.3 04/02/25 SMI"

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
 * Table of MM # constant pi/2, used in atan and k_mi_rem_pio2.
 */

#include "mcr.h"

/*
 * pio2 of nw=52, 294 Hex digits (354 decimal)
 */
const int __libmcr_TBL_pio2[] = {
	0x000001, 0x000000, 0x000001, 0x921FB5, 0x4442D1, 0x846989,
	0x8CC517, 0x01B839, 0xA25204, 0x9C1114, 0xCF98E8, 0x04177D,
	0x4C7627, 0x3644A2, 0x9410F3, 0x1C6809, 0xBBDF2A, 0x33679A,
	0x748636, 0x605614, 0xDBE4BE, 0x286E9F, 0xC26ADA, 0xDAA384,
	0x8BC90B, 0x6AECC4, 0xBCFD8D, 0xE89885, 0xD34C6F, 0xDAD617,
	0xFEB96D, 0xE80D6F, 0xDBDC70, 0xD7F6B5, 0x133F4B, 0x5D3E48,
	0x22F896, 0x3FCC92, 0x50CCA3, 0xD9C8B6, 0x7B8400, 0xF97142,
	0xC77E0B, 0x31B490, 0x6C38AB, 0xA734D2, 0x2C7F51, 0xFA499E,
	0xBF06CA, 0xBA47B9, 0x475B2C, 0x38C5E6,
};

/*
 * pio4 of nw=49
 */
const int __libmcr_TBL_pio4[] = {
	1, -1,
	0xC90FDA, 0xA22168, 0xC234C4, 0xC6628B, 0x80DC1C,
	0xD12902, 0x4E088A, 0x67CC74, 0x020BBE, 0xA63B13, 0x9B2251, 0x4A0879,
	0x8E3404, 0xDDEF95, 0x19B3CD, 0x3A431B, 0x302B0A, 0x6DF25F, 0x14374F,
	0xE1356D, 0x6D51C2, 0x45E485, 0xB57662, 0x5E7EC6, 0xF44C42, 0xE9A637,
	0xED6B0B, 0xFF5CB6, 0xF406B7, 0xEDEE38, 0x6BFB5A, 0x899FA5, 0xAE9F24,
	0x117C4B, 0x1FE649, 0x286651, 0xECE45B, 0x3DC200, 0x7CB8A1, 0x63BF05,
	0x98DA48, 0x361C55, 0xD39A69, 0x163FA8, 0xFD24CF, 0x5F8365, 0x5D23DC,
};

/*
 * Table of constants for 2/pi, used in __rem_pio2 (trig) function.
 */

/*
 * 396 Hex digits (476 decimal), 66 24-bit words of 2/pi = 0.637...
 */
const int __libmcr_TBL_ipio2[] = {
	0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62,
	0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
	0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
	0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41,
	0x3991D6, 0x398353, 0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8,
	0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
	0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
	0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08,
	0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
	0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880,
	0x4D7327, 0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B,
};
