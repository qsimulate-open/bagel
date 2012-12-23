#!/bin/sh

#g++ -std=c++11 -lblas -llapack -lgmp -lmpfr *.cc *.c ../_breitroot_3.f
g++ -std=c++11 -DBREIT -lblas -llapack -lgmp -lmpfr *.cc *.c ../_breitroot_3.f
