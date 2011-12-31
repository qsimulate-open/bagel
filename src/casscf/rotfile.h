//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_CASSCF_ROTFILE_H
#define __NEWINT_CASSCF_ROTFILE_H

#include <list>
#include <memory>
#include <algorithm>
#include <iostream>
#include <src/util/f77.h>

class RotFile {
  protected:
    double* data_;
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    int size_;
    const bool superci_;

  public:
    RotFile(const int iclos, const int iact, const int ivirt, const bool superci = true)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), superci_(superci) {
      size_ = iclos*iact+iclos*ivirt+iact*ivirt+(superci ? 1 : 0);
      data_ = new double[size_];
    };
    RotFile(const RotFile& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), size_(o.size_), superci_(o.superci_) {
      data_ = new double[size_];
      std::copy(o.data_, o.data_+size_, data_);
    };
    RotFile(const std::shared_ptr<RotFile> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), size_(o->size_), superci_(o->superci_) {
      data_ = new double[size_];
      std::copy(o->data_, o->data_+size_, data_);
    };

    ~RotFile() { delete[] data_; };

    // size of the file
    const int size() const { return size_; };
    // zero out
    void zero() { std::fill(data_, data_+size_, 0.0); };
    // returns dot product
    double ddot(RotFile& o) { return ddot_(size_, data_, 1, o.data_, 1); };
    // returns norm of the vector
    double norm() { return std::sqrt(ddot(*this)); };
    // daxpy added to self
    void daxpy(double a, RotFile& o) { daxpy_(size_, a, o.data_, 1, data_, 1); }; 
    // orthogonalize to the liset of RotFile's
    double orthog(std::list<std::shared_ptr<RotFile> > c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter)
        this->daxpy(- this->ddot(**iter), **iter);
      const double scal = 1.0/this->norm();
      dscal_(size_, scal, data_, 1);
      return 1.0/scal;
    };

    // return data_
    double* data() { return data_; };
    // return data_
    double* begin() { return data_; };
    // return data_
    double* end() { return data_+size_; };

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


    void print() {
      if (nact_ && nclosed_) {
      std::cout << " printing closed-active block" << std::endl;
      for (int i = 0; i != nact_; ++i) {
        for (int j = 0; j != nclosed_; ++j) {
          std::cout << std::setw(10) << std::setprecision(6) << ele_ca(j,i); 
        }
        std::cout << std::endl;
      }}
      if (nact_ && nvirt_) {
      std::cout << " printing virtual-active block" << std::endl;
      for (int i = 0; i != nact_; ++i) {
        for (int j = 0; j != nvirt_; ++j) {
          std::cout << std::setw(10) << std::setprecision(6) << ele_va(j,i); 
        }
        std::cout << std::endl;
      }}
      if (nclosed_ && nvirt_) {
      std::cout << " printing virtual-closed block" << std::endl;
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_; ++j) {
          std::cout << std::setw(10) << std::setprecision(6) << ele_vc(j,i); 
        }
        std::cout << std::endl;
      }}
    }; 
};


class QFile {
  protected:
    double* data_;
    int na_;
    int nb_;

  public:
    QFile(const int nb, const int na) : na_(na), nb_(nb) {
      data_ = new double[na*nb];
      std::fill(data_, data_+na*nb, 0.0);
    };
    QFile(const QFile& o) : na_(o.na_), nb_(o.nb_) {
      data_ = new double[na_*nb_];
      std::copy(o.data_, o.data_+na_*nb_, data_);
    };
    ~QFile() { delete[] data_; };

    double* data() { return data_; };
    double& element(const int i, const int j) { return data_[i+j*nb_]; };
    double* element_ptr(const int i, const int j) { return data_+i+j*nb_; };

    void zero() { std::fill(data_, data_+na_*nb_, 0.0); };

    QFile operator*(const QFile& o) {
      const int n = nb_;
      const int nc = na_; assert(nc == o.nb_);
      const int m = o.na_;
      QFile out(m, n);
      dgemm_("N", "N", n, m, nc, 1.0, data_, n, o.data_, nc, 0.0, out.data_, n); 
      return out;
    };

    QFile& operator*=(const double a) { dscal_(nb_*na_, a, data_, 1); };

    QFile operator*(const double a) {
      QFile tmp(*this);
      dscal_(nb_*na_, a, tmp.data_, 1);
      return tmp;
    };

    QFile operator/(const double a) {
      QFile tmp(*this);
      dscal_(nb_*na_, 1.0/a, tmp.data_, 1);
      return tmp;
    };

    QFile operator%(const QFile& o) {
      const int n = na_;
      const int nc = nb_; assert(nc == o.nb_);
      const int m = o.na_;
      QFile out(m, n);
      dgemm_("T", "N", n, m, nc, 1.0, data_, n, o.data_, nc, 0.0, out.data_, n); 
      return out;
    };

    QFile operator^(const QFile& o) {
      const int n = nb_;
      const int nc = na_; assert(nc == o.na_);
      const int m = o.nb_;
      QFile out(m, n);
      dgemm_("N", "T", n, m, nc, 1.0, data_, n, o.data_, nc, 0.0, out.data_, n); 
      return out;
    };


    QFile& operator+=(const QFile& o) {
      daxpy_(nb_*na_, 1.0, o.data_, 1, data_, 1); 
    };

    QFile operator+(const QFile& o) {
      QFile tmp(*this);
      daxpy_(nb_*na_, 1.0, o.data_, 1, tmp.data_, 1); 
      return tmp;
    };

    QFile operator-(const QFile& o) {
      QFile tmp(*this);
      daxpy_(nb_*na_, -1.0, o.data_, 1, tmp.data_, 1); 
      return tmp;
    };


    QFile& operator=(const QFile& o) {
      dcopy_(nb_*na_, o.data_, 1, data_, 1);
    };

    void print() const {
      std::cout << " -- printing QFile content --" << std::endl;
      const double* d = data_;
      for (int i = 0; i != na_; ++i) {
        for (int j = 0; j != nb_; ++j, ++d)
          std::cout << std::setw(12) << std::setprecision(8) << *d; 
        std::cout << std::fixed << std::endl;
      }
    };
};

#endif
