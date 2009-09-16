#ifndef __src_rysint_f77_h
#define __src_rysint_f77_h

#include <complex>

extern "C" {
 void mytranspose_complex_(const std::complex<double>*, const int*, const int*, std::complex<double>*);

};

#endif
