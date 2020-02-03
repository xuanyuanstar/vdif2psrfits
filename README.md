# vdif2psrfitsALMA

https://github.com/xuanyuanstar/vdif2psrfits

A software tool to convert ALMA VLBI data in VDIF format into pulsar search mode data in PSRFITS format.

# Installation
 * Install cfitsio and fftw library. The code has been used with cfitsio3410 and fftw3.3.4.;
 * Modify "Makefile.am" to specify route to the source code (where it finds the PSRFITs header template);
 * In the "vdif2psrfits/" directory, run "./bootstrap";
 * Run "./configure". You may need to specify links to cfitsio and fftw library manually;
 * Run "make && make install".

# Usage
 * -f   Observing central frequency (MHz)
 * -i   Input vdif pol0
 * -j   Input vdif pol1
 * -n   Band sense (-1 for lower-side band, 1 for upper-side band)
 * -s   Seconds to get statistics to fill in invalid frames
 * -k   Seconds to skip from the beginning (default 0)
 * -t   Time sample scrunch factor (by default 1). One time sample 8 microsecond
 * -S   Name of the source (by default J0000+0000)
 * -r   RA of the source
 * -c   Dec of the source
 * -D   Ouput data status (I for Stokes I, C for coherence product, X for pol0 I, Y for pol1 I, by default C)
 * -O   Route of the output file(s)
 * -h   Available options

Example of a command:
 * vdif2psrfitsALMA -f 86268.0 -i X.vdif -j Y.vdif -n -1 -O YOUR_PSRFITS/
