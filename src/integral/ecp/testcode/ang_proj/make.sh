#!/bin/sh

g++ -O3 -std=c++11 -fopenmp -lblas -llapack -lgsl -I/usr/local/include -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL main.cc -o run
#g++ -O3 -std=c++11 -fopenmp -lblas -llapack -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL ../../carsphlist.cc ../../_carsph* main.cc -o run
#g++ -Wall -Werror -O3 -std=c++11 -fopenmp -lblas -llapack -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL ../../carsphlist.cc ../../_carsph* main.cc -o run
chmod 700 run
