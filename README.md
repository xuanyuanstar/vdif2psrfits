# vdif2psrfits

https://github.com/xuanyuanstar/vdif2psrfits

Software routines to convert VLBI VDIF format into pulsar search mode data in PSRFITS format.

# Installation
 * Install cfitsio and fftw library. The code has been used with cfitsio3410 and fftw3.3.4.;
 * Modify "Makefile.am" to specify route to the source code (where it finds the PSRFITs header template);
 * In the "vdif2psrfits/" directory, run "./bootstrap";
 * Run "./configure". You may need to specify links to cfitsio and fftw library manually;
 * Run "make && make install".

# vditf2psrfitsALMA
# Dealing with ALMA VDIF output, 32 x 62.5 MHz channels
# Usage
  -f      Observing central frequency (MHz)
  -i      Input vdif pol0
  -j      Input vdif pol1
  -n      Band sense (-1 for lower-side, 1 for upper-side, by default 1)
  -s      Seconds to get statistics to fill in invalid frames
  -k      Seconds to skip from the beginning (default 0)
  -t      Time sample scrunch factor (by default 1). One time sample 8 microsecond
  -S      Name of the source (by default J0000+0000)
  -r      RA of the source (AA:BB:CC.DD)
  -c      Dec of the source (+AA:BB:CC.DD)
  -D      Ouput data status (I for Stokes I, C for coherence product, X for pol0 I, Y for pol1 I, S for Stokes, P for polarised signal, S for stokes, by default C)
  -O      Route of the output file(s).
  -h      Available options

Patching power dip options:
  -P               Replace power dip raw samples with random noise 
  -M               Replace power dip detections with mean 
  -p               Starting phase of data in scan+dip cycle (by default 0)

Example of a command:
 * vdif2psrfitsALMA -f 86268.0 -i X.vdif -j Y.vdif -n -1 -O YOUR_PSRFITS/

# vdif2psrfitsPico
# Deal with Pico & LMT VDIF output, 1 x 2 GHz channel
# Usage
  -f   Observing central frequency (MHz)
  -i   Input vdif pol0
  -j   Input vdif pol1
  -b   Band sense (-1 for lower-side, 1 for upper-side, by default 1)
  -s   Seconds to get statistics to fill in invalid frames
  -t   Time sample scrunch factor (by default 1). One time sample 8 microsecond
  -S   Name of the source (by default J0835-4510)
  -r   RA of the source (AA:BB:CC.DD)
  -c   Dec of the source (+AA:BB:CC.DD)
  -D   Ouput data status (I for Stokes I, C for coherence product, X for pol0 I, Y for pol1 I, S for Stokes, P for polarised signal, S for stokes, by default C)
  -n   Number of channels kept (Power of 2 up to 4096, by default 1)
  -O   Route of the output file 
  -h   Available options
