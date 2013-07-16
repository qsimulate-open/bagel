#!/bin/sh

#g++ -O0 -c ../_breitroot_*.cc
#g++ -O0 -c ../_spin2root_*.cc

#g++ -std=c++11 -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o -o gen
#g++ -Ofast -std=c++11 -DBREIT -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen
#g++ -Ofast -std=c++11 -DSPIN2 -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen

g++ -Ofast -std=c++11 -lblas -llapack -lgmp -lmpfr *.cc -o gen
#g++ -Ofast -std=c++11 -DBREIT -lblas -llapack -lgmp -lmpfr *.cc -o gen
#g++ -Ofast -std=c++11 -DSPIN2 -lblas -llapack -lgmp -lmpfr *.cc -o gen
chmod 700 gen
