#ifndef __src_slater_f77_h
#define __src_slater_f77_h

extern "C" {
 void root1_(const double*, const double*, double*, double*, const int*);
 void root2_(const double*, const double*, double*, double*, const int*);
 void root3_(const double*, const double*, double*, double*, const int*);
 void root4_(const double*, const double*, double*, double*, const int*);
 void root5_(const double*, const double*, double*, double*, const int*);
 void root6_(const double*, const double*, double*, double*, const int*);
 void root7_(const double*, const double*, double*, double*, const int*);
 void root8_(const double*, const double*, double*, double*, const int*);
 void root9_(const double*, const double*, double*, double*, const int*);

 void mytranspose_(const double*, const int*, const int*, double*);
 void mytranspose1_(const double*, const int*, const int*, double*);
 void mytranspose4_(const double*, const int*, const int*, double*);

 void dcopy_(const int*, const double*, const int*, double*, const int*);
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);

};

#endif
