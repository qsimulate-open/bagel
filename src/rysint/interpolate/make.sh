#!/bin/sh

g++ -lblas -llapack -lgmp -lmpfr *.cc *.c
#g++ -DBREIT -lblas -llapack -lgmp -lmpfr *.cc *.c
