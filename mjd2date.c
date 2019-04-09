// Convert from MJD to YEAR-MONTH-DAY-HOUR:MINUTE:SECOND
void mjd2date (double mjd,
			   char *date)
{
  double f, jd, dday;
  int z, alpha, a, b, c, d, e;
  int year, month, day, hour, min;
  float sec, x;

  jd = mjd + 2400000.5;
  jd += 0.5;

  z = floor (jd);
  f = fmod (jd, 1.0);

  if (z < 2299161)
    a = z;
  else
	{
	  alpha = floor ((z - 1867216.25) / 36524.25);
	  a = z + 1 + alpha - floor (alpha / 4.0);
	}
  b = a + 1524;
  c = floor ((b - 122.1) / 365.25);
  d = floor (365.25 * c);
  e = floor ((b - d) / 30.6001);

  dday = b - d - floor (30.6001 * e) + f;
  if (e < 14)
    month = e - 1;
  else
    month = e - 13;

  if (month > 2)
    year = c - 4716;
  else
    year = c - 4715;

  day = (int) floor (dday);
  x = 24.0 * (dday - day);
  x = 3600.0 * fabs (x);
  sec = fmod (x, 60.0);
  x = (x - sec) / 60.0;
  min = fmod (x, 60.0);
  x = (x - min) / 60.0;
  hour = x;

  sprintf (date, "%04d-%02d-%02d-%02d:%02d:%02.0f", year, month, day, hour, min, sec);

  return;
}
