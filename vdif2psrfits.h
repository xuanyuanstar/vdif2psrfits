void mjd2date(double mjd,char *date);
static void convert2to8(unsigned char *dest, const unsigned char *src, int bytes);
void getVDIFframetotal(const unsigned char *src, int fbytes, int *tot);
long double get_vdif_mjd(int mon, long sec);
