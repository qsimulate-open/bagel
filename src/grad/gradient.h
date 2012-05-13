//
// Author : Toru Shiozaki
// Date   : May 2012
//

#ifndef __SRC_GRAD_GRADIENT_H
#define __SRC_GRAD_GRADIENT_H

// a class for using the BFGS solver, which requires
// "clone, ddot and daxpy, along with overloaded operators and a copy constructor"

#include <vector>
#include <memory>
#include <cassert>
#include <src/util/f77.h>

class Gradient {
  protected:
    std::vector<double> data_;

  public:
    Gradient(const size_t i, const double a = 0.0) : data_(i, a) {};
    Gradient(std::vector<double> a) : data_(a.begin(), a.end()) { assert(data_.size()%3 == 0); };
    Gradient(const Gradient& o) : data_(o.data_.begin(), o.data_.end()) {};
    ~Gradient() {};

    double ddot(const Gradient& o) const { return ddot_(size(), data(), 1, o.data(), 1); };
    double ddot(const std::shared_ptr<Gradient> o) const { return ddot(*o); };

    double* data() { return &data_[0]; };
    const double* data() const { return &data_[0]; };
    size_t size() const { return data_.size(); }; 

    void daxpy(const double a, const Gradient& o) { daxpy_(size(), a, o.data(), 1, data(), 1); }; 
    void daxpy(const double a, const std::shared_ptr<Gradient> o) { daxpy(a, *o); };

    Gradient operator-(const Gradient& o) {
      Gradient out(*this);
      out.daxpy(-1.0, o); 
      return out;
    };

    double data(int i, int j) { return data_.at(3*i+j); };

};

#endif
