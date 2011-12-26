//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_CASSCF_ROTFILE_H
#define __NEWINT_CASSCF_ROTFILE_H

#include <list>
#include <memory>
#include <algorithm>
#include <src/casscf/f77.h>

class RotFile {
  protected:
    double* data_;
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    int size_;
    const bool superci_;

    // just for convenience
    static const int unit_ = 1;
  public:
    RotFile(const int iclos, const int iact, const int ivirt, const bool superci = true)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), superci_(superci) {
      size_ = iclos*iact+iclos*ivirt+iact*ivirt+(superci ? 1 : 0);
      data_ = new double[size_];
    };
    RotFile(const RotFile& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), size_(o.size_), superci_(o.superci_) {
      data_ = new double[size_];
      std::copy(o.data_, o.data_+size_, data_);
    }

    ~RotFile() { delete[] data_; };

    // size of the file
    const int size() const { return size_; };
    // zero out
    void zero() { std::fill(data_, data_+size_, 0.0); };
    // returns dot product
    double ddot(RotFile& o) { return ddot_(&size_, data_, &unit_, o.data_, &unit_); };
    // returns norm of the vector
    double norm() { return std::sqrt(ddot(*this)); };
    // daxpy added to self
    void daxpy(double a, RotFile& o) { daxpy_(&size_, &a, o.data_, &unit_, data_, &unit_); }; 
    // orthogonalize to the liset of RotFile's
    double orthog(std::list<std::shared_ptr<RotFile> > c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter)
        this->daxpy(- this->ddot(**iter), **iter);
      const double scal = 1.0/this->norm();
      dscal_(&size_, &scal, data_, &unit_);
      return 1.0/scal; 
    }

    // closed-active block. closed runs first
    double* ptr_ca() { return data_; };
    double& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; };
    // active-virtual block. virtual runs first
    double* ptr_va() { return data_ + nclosed_*nact_; };
    double& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; };
    // closed-virtual block. virtual runs first
    double* ptr_vc() { return data_ + (nclosed_+nvirt_)*nact_; };
    double& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; };
    // reference config.
    double& ele_ref() { assert(superci_); return data_[size_-1]; };

};

#endif
