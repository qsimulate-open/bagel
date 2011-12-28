//
// Author : Toru Shiozaki
// Date   : June 2009
//

#include <src/slater/srootlist.h>

extern "C" {
 // slater integrals
 void root1_(const double*, const double*, double*, double*, const int*);
 void root2_(const double*, const double*, double*, double*, const int*);
 void root3_(const double*, const double*, double*, double*, const int*);
 void root4_(const double*, const double*, double*, double*, const int*);
 void root5_(const double*, const double*, double*, double*, const int*);
 void root6_(const double*, const double*, double*, double*, const int*);
 void root7_(const double*, const double*, double*, double*, const int*);
 void root8_(const double*, const double*, double*, double*, const int*);
 void root9_(const double*, const double*, double*, double*, const int*);
}

SRootList::SRootList() {

  srfunc[1] = &root1_;
  srfunc[2] = &root2_;
  srfunc[3] = &root3_;
  srfunc[4] = &root4_;
  srfunc[5] = &root5_;
  srfunc[6] = &root6_;
  srfunc[7] = &root7_;
  srfunc[8] = &root8_;
  srfunc[9] = &root9_;
}


SRootList::~SRootList() {

}
