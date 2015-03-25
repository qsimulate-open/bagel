#!/bin/sh

g++ -std=c++11 -lblas -llapack -lgmp -lmpfr dfact.cc -o dfact
chmod 700 dfact
