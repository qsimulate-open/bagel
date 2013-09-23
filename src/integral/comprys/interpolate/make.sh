#!/bin/sh

g++ -O3 -std=c++11 -lblas -llapack -L/opt/local/lib -lgmp -L/../../ -lmpfr *.cc -o gen
chmod 700 gen
