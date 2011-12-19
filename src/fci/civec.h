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
    std::shared_ptr<std::vector<double> > cc_;
    int lena_;
    int lenb_;

  public:
    Civec(const size_t lb, const size_t la) : lena_(la), lenb_(lb) {
      std::shared_ptr<std::vector<double> > tmp(new std::vector<double>); 
      tmp->resize(la*lb);
      std::fill(tmp->begin(), tmp->end(), 0.0);
      cc_ = tmp;
    }
    ~Civec() {};

    double& element(size_t i, size_t j) { return (*cc_)[i+j*lenb_]; }; // I RUNS FIRST 
    double* element_ptr(size_t i, size_t j) { return &((*cc_)[i+j*lenb_]); }; // I RUNS FIRST 

    void zero() { std::fill(cc_->begin(), cc_->end(), 0.0); };

    std::shared_ptr<Civec> transpose() {
      std::shared_ptr<Civec> ct(new Civec(lena_, lenb_));
      std::shared_ptr<std::vector<double> > cct = ct->cc(); 
      mytranspose_(&((*cc_)[0]), &lenb_, &lena_, &((*cct)[0])); 
      return ct;
    };

    std::shared_ptr<std::vector<double> > cc() { return cc_; };
};

#endif
