//
// Author : Toru Shiozaki
// Date   : June 2009
//

#ifndef __src_rysint_erirootlist_h
#define __src_rysint_erirootlist_h

#include <src/rysint/macros.h>

struct ERIRootList  {
  public:
    ERIRootList();
    ~ERIRootList();

    void (*rfunc[RYS_MAX + 1])(const double*, double*, double*, const int*);

    void rootfunc_call(const unsigned int i, const double* a1, double* a2, double* a3, const int* a4) {
      return (rfunc[i])(a1, a2, a3, a4);
    };

}; 

#endif

