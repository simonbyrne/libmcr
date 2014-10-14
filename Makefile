#
# Copyright 04/09/04 Sun Microsystems, Inc.  All Rights Reserved.
#

# libmcr -- libm correctly rounded
#
#

# !!!!!!
# Installation specific
##############################################
CC		= cc
OPT		= -g
CFLAGS		=
##############################################

SRCDIR		= ./src
TESTSDIR	= ./tests


all: src tests

runtests: FRC
	(cd ${TESTSDIR}; make all CC=${CC} OPT=${OPT})

src: FRC
	(cd ${SRCDIR}; make CC=${CC} OPT=${OPT} CFLAGS=${CFLAGS})

tests: FRC
	(cd ${TESTSDIR}; make CC=${CC} OPT=${OPT})

clean:
	(cd ${SRCDIR}; make clean)
	(cd ${TESTSDIR}; make clean)

clobber:
	(cd ${SRCDIR}; make clobber)
	(cd ${TESTSDIR}; make clobber)

FRC:

