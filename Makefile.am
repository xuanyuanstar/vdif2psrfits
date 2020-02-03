AM_CFLAGS = -I@top_srcdir@/ @CFITSIO_CFLAGS@ @FFTW_CFLAGS@ 

#bin_PROGRAMS= vdif2psrfitsPico vdif2psrfitsALMA
bin_PROGRAMS = vdif2psrfitsALMA vdif2psrfitsPico
lib_LTLIBRARIES = libVDIF.la

libVDIF_la_SOURCES = dec2hms.c polyco.c vdifio.c write_psrfits.c cvrt2to8.c mjd2date.c getVDIFFrameDetection.c downsample.c
libVDIF_la_LIBADD = @CFITSIO_LIBS@ @FFTW_LIBS@ 

vdif2psrfitsALMA_SOURCES = vdif2psrfitsALMA.c
vdif2psrfitsALMA_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@

vdif2psrfitsPico_SOURCES = vdif2psrfitsPico.c
vdif2psrfitsPico_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@

AM_CPPFLAGS = -DPSRFITS_TEMPLATE_DIR='"/cluster/pulsar/kliu/Soft/vdif2psrfits/"'
