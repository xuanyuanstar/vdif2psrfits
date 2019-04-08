// Convert 2-bit string to 8-bit, ranging from 0 to 3
// Arguments are output, input, number of bytes in input
// Routine adopted from vdifio

static void convert2to8(unsigned char *dest, const unsigned char *src, int bytes)
{
  const unsigned char levels[4] = {0, 1, 2, 3};
  static int first = 1;
  static unsigned char lut2to8[256][4];
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
