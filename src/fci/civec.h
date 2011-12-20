//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef NEWINT_FCI_CIVEC_H
#define NEWINT_FCI_CIVEC_H

#include <memory>
#include <vector>
#include <algorithm>
#include <src/fci/f77.h>

class Civec {
  protected:
    // !!CAUTION!!
    // cc is formated so that B runs first.
    double* cc_;
    int lena_;
    int lenb_;
    const bool alloc;

  public:
    Civec(const size_t lb, const size_t la) : lena_(la), lenb_(lb), alloc(true) {
      cc_ = new double[la*lb];
      std::fill(cc_, cc_ + la*lb, 0.0);
    }
    Civec(const size_t lb, const size_t la, double* din_) : lena_(la), lenb_(lb), alloc(false) {
      cc_ = din_; 
      std::fill(cc_, cc_ + la*lb, 0.0);
    }
    ~Civec() {
      if (alloc) delete[] cc_;
    };

    double& element(size_t i, size_t j) { return cc_[i+j*lenb_]; }; // I RUNS FIRST 
    double* element_ptr(size_t i, size_t j) { return cc_+i+j*lenb_; }; // I RUNS FIRST 
    double* first() { return cc_; };

    void zero() { std::fill(cc_, cc_+lena_*lenb_, 0.0); };

    std::shared_ptr<Civec> transpose() {
      std::shared_ptr<Civec> ct(new Civec(lena_, lenb_));
      double* cct = ct->first(); 
      mytranspose_(cc_, &lenb_, &lena_, cct); 
      return ct;
    };
    int lena() const { return lena_; };
    int lenb() const { return lenb_; };

    // some functions for convenience
    double ddot(Civec& other) {
      const int lab = lena_ * lenb_;
      const int unit = 1;
      return ddot_(&lab, cc_, &unit, other.first(), &unit);
    };
    void daxpy(double a, Civec& other) {
      const int lab = lena_ * lenb_;
      const int unit = 1;
      daxpy_(&lab, &a, other.first(), &unit, cc_, &unit);
    };
    double norm() {
      const int lab = lena_ * lenb_;
      const int unit = 1;
      return std::sqrt(ddot_(&lab, cc_, &unit, cc_, &unit));
    };
};

class Dvec {
  protected:
    std::vector<std::shared_ptr<Civec> > dvec_;
    double* data_;
    int lenb_;
    int lena_;
    int ij_;
  public:
    Dvec(const size_t lb, const size_t la, const size_t ij) : lena_(la), lenb_(lb), ij_(ij) {
      // actually data should be in a consecutive area to call dgemm.
      data_ = new double[lb*la*ij];
      double* tmp = data_;
      for (int i = 0; i != ij; ++i, tmp+=lb*la) {
        std::shared_ptr<Civec> c(new Civec(lb, la, tmp)); 
        dvec_.push_back(c);
      }
    };
    ~Dvec() {
      delete[] data_;
    };
    std::shared_ptr<Civec> data(const size_t i) { return dvec_.at(i); };
    void zero() { std::fill(data_, data_+lena_*lenb_*ij_, 0.0); };
    double* first() { return data_; };

    int lena() const { return lena_; };
    int lenb() const { return lenb_; };
    int ij() const { return ij_; };
};

#endif
