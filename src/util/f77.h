#ifndef __src_rysint_f77_h
#define __src_rysint_f77_h

#include <complex>

extern "C" {
 void mytranspose_complex_(const std::complex<double>*, const int*, const int*, std::complex<double>*);
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*); 
 const double ddot_(const int*, const double*, const int*, const double*, const int*); 

};

#endif
