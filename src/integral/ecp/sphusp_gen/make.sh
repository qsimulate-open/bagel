#!/bin/sh

g++ -Wall -Werror -O3 -std=c++11 -fopenmp -lblas -llapack -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL ../sphharmonics.cc sphusp_gen.cc -o gen
chmod 700 gen
