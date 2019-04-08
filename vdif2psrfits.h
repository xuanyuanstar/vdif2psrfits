#include "mjd2date.c"
#include "cvrt2to8.c"

void mjd2date(double mjd,char *date);
static void convert2to8(unsigned char *dest, const unsigned char *src, int bytes);
