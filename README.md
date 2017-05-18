#Start with config/, put cfitsio.m4 and a few other files in it

#Create configure based on configure.ac, grap more files into config/

./bootstrap

#Create Makefile based on Makefile.am

./configure

#Create output

make
