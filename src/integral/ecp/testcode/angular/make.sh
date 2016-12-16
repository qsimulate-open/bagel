#!/bin/sh

g++ -O0 -std=c++11 -fopenmp -lblas -llapack -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL ../../../../util/math/sphharmonics.cc ../../_sphusp_*.cc main.cc -o run
chmod 700 run
