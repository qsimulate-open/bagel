//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_FCI_RDM_H
#define __NEWINT_FCI_RDM_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <src/util/f77.h>

template <int rank>
class RDM {
  protected:
    double* data_;
    const int norb_;
    int dim_;

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
    double* data() { return data_; };
    double& element(int i, int j) { return *(data_+i+j*dim_); };
    const double& element(int i, int j) const { return *(data_+i+j*dim_); };
    const double* element_ptr(int i, int j) const { return data_+i+j*dim_; };
    // careful, this should not be called for those except for 2RDM. 
    double& element(int i, int j, int k, int l) { assert(rank == 2); return *(data_+i+norb_*(j+norb_*(k+norb_*l))); };
    double* element_ptr(int i, int j, int k, int l) { assert(rank == 2); return data_+i+norb_*(j+norb_*(k+norb_*l)); };
    const double& element(int i, int j, int k, int l) const { assert(rank == 2); return *(data_+i+norb_*(j+norb_*(k+norb_*l))); };

    void zero() { std::fill(data_, data_+dim_*dim_, 0.0); };
    void daxpy(const double a, const RDM& o) {
      const int size = dim_*dim_;
      const int unit = 1;
      daxpy_(&size, &a, o.data_, &unit, data_, &unit);
    };

    void daxpy(const double a, const std::shared_ptr<RDM>& o) { this->daxpy(a, *o); };

    std::vector<double> diag() const {
      std::vector<double> out(dim_);
      for (int i = 0; i != dim_; ++i) out[i] = element(i,i);
      return out;
    };

    std::pair<std::vector<double>, std::vector<double> > generate_natural_orbitals() const {
      assert(rank == 1);
      std::vector<double> buf(dim_*dim_);
      std::vector<double> vec(dim_);
#define ALIGN
#ifdef ALIGN
      for (int i = 0; i != dim_; ++i) buf[i+i*dim_] = 2.0; 
      daxpy_(dim_*dim_, -1.0, data_, 1, &(buf[0]), 1);
#else
      daxpy_(dim_*dim_, 1.0, data_, 1, &(buf[0]), 1);
#endif
      int lwork = 5*dim_;
      std::vector<double> work(lwork);
      int info;
      dsyev_("V", "U", &dim_, &(buf[0]), &dim_, &(vec[0]), &(work[0]), &lwork, &info);
      assert(!info);
#ifdef ALIGN
      for (auto i = vec.begin(); i != vec.end(); ++i) *i = 2.0-*i;
#endif
      return std::make_pair(buf, vec);
    };

    void transform(const std::vector<double>& coeff) {
      const double* start = &(coeff[0]);
      double* buf = new double[dim_*dim_];
      if (rank == 1) {
        dgemm_("N", "N", dim_, dim_, dim_, 1.0, data_, dim_, start, dim_, 0.0, buf,   dim_); 
        dgemm_("T", "N", dim_, dim_, dim_, 1.0, start, dim_, buf,   dim_, 0.0, data_, dim_); 
      } else if (rank == 2) {
        // first half transformation
        dgemm_("N", "N", dim_*norb_, norb_, norb_, 1.0, data_, dim_*norb_, start, norb_, 0.0, buf, dim_*norb_); 
        for (int i = 0; i != norb_; ++i)
          dgemm_("N", "N", dim_, norb_, norb_, 1.0, buf+i*dim_*norb_, dim_, start, norb_, 0.0, data_+i*dim_*norb_, dim_); 
        // then tranpose
        mytranspose_(data_, &dim_, &dim_, buf);
        // and do it again
        dgemm_("N", "N", dim_*norb_, norb_, norb_, 1.0, buf, dim_*norb_, start, norb_, 0.0, data_, dim_*norb_); 
        for (int i = 0; i != norb_; ++i)
          dgemm_("N", "N", dim_, norb_, norb_, 1.0, data_+i*dim_*norb_, dim_, start, norb_, 0.0, buf+i*dim_*norb_, dim_); 
        // to make sure for non-symmetric density matrices (and anyway this should be cheap).
        mytranspose_(buf, &dim_, &dim_, data_);
      } else {
        assert(false);
      }
      delete[] buf;
    };

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
                if (std::abs(element(l,k,j,i)) > 1.0e-0) std::cout << std::setw(3) << l << std::setw(3)
                      << k << std::setw(3) << j << std::setw(3) << i
                      << std::setw(12) << std::setprecision(7) << element(l,k,j,i) << std::endl;
        } } } }
      }
    };
};

#endif
