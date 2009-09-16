#ifndef __src_rysint_f77_h
#define __src_rysint_f77_h

extern "C" {
 void rysroot_(const double*, double*, double*, const int*, const int*);
 void mytranspose_(const double*, const int*, const int*, double*);
 void mytranspose1_(const double*, const int*, const int*, double*);
 void mytranspose4_(const double*, const int*, const int*, double*);

 void dcopy_(const int*, const double*, const int*, double*, const int*);
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);

 void eriroot0_(const double*, double*, double*, const int*);
 void eriroot1_(const double*, double*, double*, const int*);
 void eriroot2_(const double*, double*, double*, const int*);
 void eriroot3_(const double*, double*, double*, const int*);
 void eriroot4_(const double*, double*, double*, const int*);
 void eriroot5_(const double*, double*, double*, const int*);
 void eriroot6_(const double*, double*, double*, const int*);
 void eriroot7_(const double*, double*, double*, const int*);
 void eriroot8_(const double*, double*, double*, const int*);
 void eriroot9_(const double*, double*, double*, const int*);
};

#endif
