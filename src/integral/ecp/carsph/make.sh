#!/bin/sh

g++ -Wall -Werror -O3 -std=c++11 -fopenmp -lblas -llapack -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL carsph_tab.cc -o tab
chmod 700 tab
