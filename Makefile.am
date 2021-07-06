ACLOCAL_AMFLAGS = -I config
AM_CFLAGS = -I@top_srcdir@/ @CFITSIO_CFLAGS@ @FFTW_CFLAGS@ -DPSRFITS_TEMPLATE_DIR='"/home/kliu/Soft/MPIvdif2psrfits/templates/"'

bin_PROGRAMS = MPIvdif2psrfitsPico MPIvdif2psrfitsALMA
lib_LTLIBRARIES = libVDIF.la

libVDIF_la_SOURCES = dec2hms.c polyco.c vdifio.c write_psrfits.c cvrt2to8.c mjd2date.c getVDIFFrameDetection.c downsample.c
libVDIF_la_LIBADD = @CFITSIO_LIBS@ @FFTW_LIBS@ 

MPIvdif2psrfitsPico_SOURCES = MPIvdif2psrfitsPico.c
MPIvdif2psrfitsPico_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f

MPIvdif2psrfitsALMA_SOURCES = MPIvdif2psrfitsALMA.c
MPIvdif2psrfitsALMA_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f
