//Convert 2-bit string to 8-bit, ranging 0-255
//Arguments are output, input, number of bytes in input
static void convert2to8(unsigned char *dest, const unsigned char *src, int bytes)
{
  /* choose levels such that the ratio of high to low is as close to 3.3359      
   * as possible to best maintain amplitude scaling.  127.5 is the center of  
   * the scale (equates to 0).  118.5/35.5 is pretty close to optimal.
   */
  const unsigned char levels[4] = {9, 92, 163, 246};
  //const unsigned char levels[4] = {0, 1, 2, 3};
  static int first = 1;
  static unsigned char lut2to8[256][4];   /* mapping from input 8 bits (4 samples) to output 4 8-bit samples */
  int i, o;

  if(first)
	{
	  /* assemble look up table */

	  for(i = 0; i < 256; ++i)
		{
		  int j;

		  for(j = 0; j < 4; ++j)
			{
			  int k;

			  k = (i >> (2*j)) & 0x3;
			  lut2to8[i][j] = levels[k];
			}
		}
	}

  o = 0;
  for(i = 0; i < bytes; ++i)
	{
	  dest[o++] = lut2to8[src[i]][0];
	  dest[o++] = lut2to8[src[i]][1];
	  dest[o++] = lut2to8[src[i]][2];
	  dest[o++] = lut2to8[src[i]][3];
	}
}
