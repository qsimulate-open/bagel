//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_FCI_RDM_H
#define __NEWINT_FCI_RDM_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/fci/f77.h>

template <int rank>
class RDM {
  protected:
    double* data_;
    const int norb_;
    size_t dim_;

  public:
    RDM(const int n) : norb_(n) {
      assert(rank > 0);
      dim_ = 1;
      for (int i = 0; i != rank; ++i) dim_ *= n;
      data_ = new double[dim_*dim_];
    };
    RDM(const RDM& o) : norb_(o.norb_), dim_(o.dim_) {
      data_ = new double[dim_*dim_];
      std::copy(o.data_, o.data_+dim_*dim_, data_);
    };
    ~RDM() { delete[] data_; };

    double* first() { return data_; };
    double& element(int i, int j) { return *(data_+i+j*dim_); };
    // careful, this should not be called for those except for 2RDM. 
    double& element(int i, int j, int k, int l) { assert(rank == 2); return *(data_+i+norb_*(j+norb_*(k+norb_*l))); };

    void zero() { std::fill(data_, data_+dim_*dim_, 0.0); };
    void daxpy(const double a, const RDM& o) {
      const int size = dim_*dim_;
      const int unit = 1;
      daxpy_(&size, &a, o.data_, &unit, data_, &unit);
    };

    void daxpy(const double a, const std::shared_ptr<RDM>& o) { this->daxpy(a, *o); };

    void print() {
      if (rank == 1) {
        for (int i = 0; i != norb_; ++i) {
          for (int j = 0; j != norb_; ++j)
            std::cout << std::setw(12) << std::setprecision(7) << element(j,i); 
          std::cout << std::endl;
        }
      } else if (rank == 2) {
        for (int i = 0; i != norb_; ++i) {
          for (int j = 0; j != norb_; ++j) {
            for (int k = 0; k != norb_; ++k) {
              for (int l = 0; l != norb_; ++l) {
                if (std::abs(element(l,k,j,i)) > 1.0e-2) std::cout << std::setw(3) << l << std::setw(3)
                      << k << std::setw(3) << j << std::setw(3) << i
                      << std::setw(12) << std::setprecision(7) << element(l,k,j,i) << std::endl;
        } } } }
      }
    };
};

#endif
