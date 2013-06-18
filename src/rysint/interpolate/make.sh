#!/bin/sh

#gfortran -c ../_breitroot_*.f
#gfortran -c ../_spin2root_*.f

#g++ -std=c++11 -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o -o gen

#g++ -Ofast -std=c++11 -DBREIT -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen
#g++ -Ofast -std=c++11 -DSPIN2 -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen
g++ -Ofast -std=c++11 -lgfortran -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen
chmod 700 gen
