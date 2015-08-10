#!/bin/sh

g++ -O3 -std=c++11 -lblas -llapack -lgsl -I/usr/local/include -L/opt/local/lib -lgmp -lmpfr -I$HOME/develop/BAGEL main.cc -o run
chmod 700 run
