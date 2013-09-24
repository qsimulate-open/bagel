#!/bin/sh

g++ -O3 -std=c++11 -I../../../../ -lblas -llapack -L/opt/local/lib -lgmp -lmpfr *.cc ../*.cc -o gen
#g++ -O3 -std=c++11 -I/Users/reynolds/develop/BAGEL/ -lblas -llapack -L/opt/local/lib -lgmp -lmpfr *.cc ../*.cc -o gen
chmod 700 gen
