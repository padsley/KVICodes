# KVICodes
KVI Codes from Muhsin Harakeh for computing e.g. form factors, energy-weighted sum rules etc.

More documentation will be coming hopefully soon.

I have to modify the codes a little bit to get them working (mainly just old FORTRAN not getting on with gfortran) and more detailed instructions may be required if you would like to use these codes.

Email me if you need help.

to compile belgen:

gfortran --free-form belgen.for bellib.o mnhlib.o cio/cio*.o -o belgen