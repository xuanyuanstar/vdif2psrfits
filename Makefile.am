ACLOCAL_AMFLAGS = -I config
AM_CFLAGS = -I@top_srcdir@/ @CFITSIO_CFLAGS@ @FFTW_CFLAGS@ -DPSRFITS_TEMPLATE_DIR='"/home/kliu/Soft/vdif2psrfits/templates/"'

bin_PROGRAMS = vdif2psrfitsALMA vdif2psrfitsPico vdif2psrfitsPicoMPI
lib_LTLIBRARIES = libVDIF.la

libVDIF_la_SOURCES = dec2hms.c polyco.c vdifio.c write_psrfits.c cvrt2to8.c mjd2date.c getVDIFFrameDetection.c downsample.c
libVDIF_la_LIBADD = @CFITSIO_LIBS@ @FFTW_LIBS@ 

vdif2psrfitsALMA_SOURCES = vdif2psrfitsALMA.c
vdif2psrfitsALMA_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f_threads

vdif2psrfitsPico_SOURCES = vdif2psrfitsPico.c
vdif2psrfitsPico_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f_threads

vdif2psrfitsPicoMPI_SOURCES = vdif2psrfitsPicoMPI.c
vdif2psrfitsPicoMPI_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f

AM_CPPFLAGS = -I@top_srcdir@/ @CFITSIO_CFLAGS@ @FFTW_CFLAGS@ -DPSRFITS_TEMPLATE_DIR='"/home/kliu/Soft/vdif2psrfits/templates/"'
