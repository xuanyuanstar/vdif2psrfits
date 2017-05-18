AC_INIT([VDIF2psrfits], [1.0], [kliu@mpifr-bonn.mpg.de])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall foreign])

AM_PROG_AR

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_PROG_F77_C_O

AC_PROG_INSTALL
AC_PROG_LIBTOOL

AC_CHECK_LIB([m],[cos])

# Checks for essential libraries
SWIN_LIB_FFTW
SWIN_LIB_CFITSIO

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 ])
 AC_OUTPUT
