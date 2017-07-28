# Autoconf requirements
# Autoconf version >= 2.69
AC_PREREQ([2.69])                      
AC_INIT([proto],[2.0.0],[urrfinjuss@gmail.com], [https://github.com/urrfinjuss/gas-network/])

PKG_PROG_PKG_CONFIG

AC_CONFIG_SRCDIR([src/comps.c],[src/init.c],[src/main.c],[src/network.c],[src/output.c],[src/verify.c],[src/header.h],[src/input.c],[src/message.c],[src/nodes.c],[src/pipes.c])
#AC_CONFIG_HEADERS([config.h])
# Initialize libtool
#LT_INIT
#Auxiliary file go here
AC_CONFIG_AUX_DIR([build-aux])
#AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_MACRO_DIRS([m4])

# Enable "automake"
AM_INIT_AUTOMAKE([1.15 subdir-objects])

# Look for a C compiler
AC_PROG_CC

# Check for libraries
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([Could not find -lm])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR([Could not find -lgslcblas])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([Could not find -lgsl])])

#Check for library functions.

# Check for header files
AC_CHECK_HEADERS([gsl/gsl_errno.h],[], [], [[#ifdef HAVE_GSL_GSL_ERRNO_H
#include <gsl/gsl_errno.h>
#endif]])
AC_CHECK_HEADERS([gsl/gsl_spline.h],[],[], [[#ifdef HAVE_GSL_GSL_SPLINE_H
#include <gsl/gsl_spline.h>
#endif]])



AC_CONFIG_FILES([Makefile])           
AC_OUTPUT
