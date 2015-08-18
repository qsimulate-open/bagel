#!/bin/sh
g++ -O0 -g -std=c++11 -lblas -llapack -lgsl -I/usr/local/include -L/opt/local/lib -I$HOME/develop/BAGEL -I. ../../../integral/rys/_eriroot_*.cc ../../../util/parallel/process.cc ../../../util/parallel/mpi_interface.cc ../../../util//f77_interface.cc *.cc -o run
chmod 700 run
