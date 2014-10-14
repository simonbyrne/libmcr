#
# Copyright 03/05/04 Sun Microsystems, Inc.  All Rights Reserved.
#

# libmcr -- libm correctly rounded
#
# Correctly Rounded elementary functions(wrap up routines)
#	double __libmcr_exp(double x)
#	double __libmcr_log(double x)
#	double __libmcr_pow(double x)
#	double __libmcr_atan(double x)
#	double __libmcr_sin(double x)
#	double __libmcr_cos(double x)
#	double __libmcr_tan(double x)
#
# Elementary functions with correction term(~20 bit extra precision)
#	double __libmcr_mx_exp(double x, double *err)
#	double __libmcr_mx_log(double x, double *err)
#	double __libmcr_mx_pow(double x, double y, double *err)
#	double __libmcr_mx_atan(double x, double *err)
#	double __libmcr_mx_sin(double x, double *err)
#	double __libmcr_mx_cos(double x, double *err)
#	double __libmcr_mx_tan(double x, double *err)
#
# Auxilary and kernel functions
#	double __libmcr_k_mx_exp(double x, double *err, int *expn)
#	double __libmcr_k_mx_pow(double *x, double *y, double *err, int *expn)
#	double __libmcr_k_mx_tan(double x, double t, double *err)
#	double __libmcr_k_mx_cos(double x, double t, double *err)
#	double __libmcr_k_mx_sin(double x, double t, double *err)
#	int    __libmcr_k_mx_rem_pio2(double x, double *y)
#	int    __libmcr_mx_check(int hz, int lz, int hr, int rnd)
#	void   __libmcr_mx_poly(double *z, double *a, double *e, int n)
#
# Multi-precision kernel functions with test flag
#	double __libmcr_mi_atan(double, int, int *, int);
#	double __libmcr_mi_exp(double, int, int *, int);
#	double __libmcr_mi_log(double, int, int *, int);
#	double __libmcr_mi_pow(double, double, int, int *, int);
#	double __libmcr_mi_sin(double, int, int *, int);
#	double __libmcr_mi_cos(double, int, int *, int);
#	double __libmcr_mi_tan(double, int, int *, int);
#
#	void __libmcr_mm_sub(int *, int *, int *, int);
#	void __libmcr_mm_add(int *, int *, int *, int);
#	void __libmcr_mm_div(int *, int *, int *, int);
#	void __libmcr_mm_divi(int *, int, int *, int);
#	void __libmcr_mm_mul(int *, int *, int *, int);
#	void __libmcr_mm_muli(int *, int, int *, int);
#	void __libmcr_mm_scalbn(int *, int, int);
#	void __libmcr_mm_sqrt(int *, int *, int);
#	void __libmcr_mm_atan(int *, int *, int);

# Auxilary and kernel functions
#	void   __libmcr_mi_dtomi(double, int *, int);
#	void   __libmcr_mi_itomi(int, int *, int);
#	void   __libmcr_mi_format(int *, int);
#	double __libmcr_mi_mitod(int *, int, int *);
#	double __libmcr_mi_final(int *, int, int *, int);
#	int    __libmcr_mi_ilogb(int *);
#	void __libmcr_mm_pio2(int *, int);
#	void __libmcr_mm_pio4(int *, int);
#	void __libmcr_mm_ln2(int *, int);
#	void __libmcr_k_mm_exp(int *, int *, int);
#	void __libmcr_k_mm_sin(int *, int *, int);
#	void __libmcr_k_mm_cos(int *, int *, int);
#	int  __libmcr_k_mi_rem_pio2(double, int *, int);
#	int  __libmcr_k_rem_pio2m(double *, double *, int, int, int);
#	void __libmcr_k_mi_log(double, int *, int);
#	double __libmcr_k_exactsqrt(double, int *);
#
# Data
#	const int __libmcr_TBL_atan_hi[];
#	const int __libmcr_TBL_atan_lo[];
#	const int __libmcr_TBL_pio2[];
#	const int __libmcr_TBL_pio4[];
#	const int __libmcr_TBL_ipio2[];
#

CC	= cc
CFLAGS	=

INCDIR	= ../include
HSRC	= ${INCDIR}/libmcr.h ${INCDIR}/mcr.h

LIBDIR	= ../lib
LIB 	= ${LIBDIR}/libmcr.a

COMPLINEOPTS = -I${INCDIR} ${OPT} ${CFLAGS}


OBJS 	= \
	__libmcr_exp.o \
	__libmcr_log.o \
	__libmcr_pow.o \
	__libmcr_atan.o \
	__libmcr_sin.o \
	__libmcr_cos.o \
	__libmcr_tan.o \
	__libmcr_mx_exp.o \
	__libmcr_mx_log.o \
	__libmcr_mx_pow.o \
	__libmcr_mx_atan.o \
	__libmcr_mx_sin.o \
	__libmcr_mx_cos.o \
	__libmcr_mx_tan.o \
	__libmcr_k_mx_exp.o \
	__libmcr_k_mx_pow.o \
	__libmcr_k_mx_tan.o \
	__libmcr_k_mx_cos.o \
	__libmcr_k_mx_sin.o \
	__libmcr_k_mx_rem_pio2.o \
	__libmcr_k_exactsqrt.o \
	__libmcr_mx_check.o \
	__libmcr_mx_poly.o \
	__libmcr_mi_atan.o \
	__libmcr_mi_exp.o \
	__libmcr_mi_log.o \
	__libmcr_mi_pow.o \
	__libmcr_mi_sin.o \
	__libmcr_mi_cos.o \
	__libmcr_mi_tan.o \
	__libmcr_mm_sub.o \
	__libmcr_mm_add.o \
	__libmcr_mm_div.o \
	__libmcr_mm_divi.o \
	__libmcr_mm_mul.o \
	__libmcr_mm_muli.o \
	__libmcr_mm_scalbn.o \
	__libmcr_mm_sqrt.o \
	__libmcr_mm_atan.o \
	__libmcr_mi_dtomi.o \
	__libmcr_mi_itomi.o \
	__libmcr_mi_format.o \
	__libmcr_mi_mitod.o \
	__libmcr_mi_final.o \
	__libmcr_mi_ilogb.o \
	__libmcr_mm_pio2.o \
	__libmcr_mm_ln2.o \
	__libmcr_k_mm_exp.o \
	__libmcr_k_mm_sin.o \
	__libmcr_k_mm_cos.o \
	__libmcr_k_mi_rem_pio2.o \
	__libmcr_k_rem_pio2m.o \
	__libmcr_k_mi_log.o \
	__libmcr_TBL_atan.o \
	__libmcr_TBL_pio2.o


all: $(LIB)

$(LIB): $(HSRC) $(OBJS)
	ar cru $(LIB) $(OBJS)

$(OBJS): Makefile

.c.o: $*.c
	$(CC) ${COMPLINEOPTS} -c -o $*.o $*.c

clean:
	/bin/rm -f $(OBJS)

clobber: clean
	rm -f ${LIB}

