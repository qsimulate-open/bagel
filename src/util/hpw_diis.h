//
// Author : Toru Shiozaki
// Date   : Jan 1st, 2012, 0.46am :-)
//

// implements special class of DIIS
// Hampel, Peterson, Werner, Chem. Phys. Lett. 190, 1 (1992).
//
// In addition to the requirement in diis.h, T should have
// functions exp(), log(), unit() that return shared_ptr<T>.
// Also operator "*=" should be overloaded in T.
// Copy constructor is needed as well.

#ifndef __NEWINT_UTIL_HPW_DIIS_H
#define __NEWINT_UTIL_HPW_DIIS_H

#include <src/util/diis.h>

template<class T>
class HPW_DIIS  {
  typedef std::shared_ptr<T> RefT;
  protected:
    DIIS<T> diis_;
    RefT base_;
    RefT prev_;
    const RefT orig_;

  public:
    HPW_DIIS(const int n, RefT o) : diis_(n), base_(o->clone()), orig_(o) { base_->unit(); };
    ~HPW_DIIS() {};

    RefT extrapolate(const RefT rot) {
      // prev = log(base)
      RefT expo = (*base_**rot).log();
      RefT prev_ = base_->log();
      RefT err(new T(prev_ ? (*expo-*prev_) : (*expo)));

      RefT extrap = diis_.extrapolate(std::make_pair(expo, err)); 
      // returns unitary matrix with respect to the original matrix
      RefT out(new T(*orig_* *extrap->exp()));
      *base_ = *extrap->exp();
      return out;
    };

};

#endif
