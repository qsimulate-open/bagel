#ifndef __src_fci_f77_h
#define __src_fci_f77_h

extern "C" {
 void mytranspose_(const double*, const int*, const int*, double*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);
};

#endif
