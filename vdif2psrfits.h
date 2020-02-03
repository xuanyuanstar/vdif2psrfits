#include "mjd2date.c"
#include "cvrt2to8.c"
double mean_2to8 = 1.5;
void mjd2date(double mjd,char *date);
static void convert2to8(unsigned char *dest, const unsigned char *src, int bytes);
void dec2hms(char *out, double in, int sflag);
