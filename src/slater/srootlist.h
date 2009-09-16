//
// Author : Toru Shiozaki
// Date   : June 2009
//

#ifndef __src_slater_srootlist_h
#define __src_slater_srootlist_h

#include <src/rysint/macros.h>

struct SRootList  {
  public:
    SRootList();
    ~SRootList();

    void (*srfunc[RYS_MAX + 1])(const double*, const double*, double*, double*, const int*);

    void srootfunc_call(const unsigned int i, const double* a0, const double* a1, double* a2, double* a3, const int* a4) {
      return (srfunc[i])(a0, a1, a2, a3, a4);
    };

}; 

#endif

