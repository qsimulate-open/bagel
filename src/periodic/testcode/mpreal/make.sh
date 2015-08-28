#!/bin/sh
g++ -O0 -g -std=c++11 -lblas -llapack -lgsl -I/usr/local/include -lgmp -lmpfr -L/opt/local/lib *.cc -o run
chmod 700 run
