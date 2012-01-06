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
#include <iomanip>
#include <cassert>
#include <src/scf/matrix1e.h>
#include <src/util/f77.h>

class RotFile {
  protected:
    std::unique_ptr<double[]> data_;
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const int size_;
    const bool superci_;

  public:
    RotFile(const int iclos, const int iact, const int ivirt, const bool superci = true)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), superci_(superci), size_(iclos*iact+iclos*ivirt+iact*ivirt+(superci ? 1 : 0)) {
        std::unique_ptr<double[]> tmp(new double[size_]);
        data_ = std::move(tmp);
    };
    RotFile(const RotFile& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), size_(o.size_), superci_(o.superci_), data_(new double[o.size_]) {
      std::cout << "size " << size_ <<  std::endl;
      std::copy(o.data(), o.data()+size_, data());
    };
    RotFile(const std::shared_ptr<RotFile> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), size_(o->size_), superci_(o->superci_), data_(new double[o->size_]) {
      std::cout << "a size " << size_ <<  std::endl;
      std::copy(o->data(), o->data()+size_, data());
    };

    ~RotFile() {  };

    // size of the file
    const int size() const { return size_; };
    // zero out
    void zero() { std::fill(data(), data()+size_, 0.0); };
    // returns dot product
    double ddot(const RotFile& o) { return ddot_(size_, data(), 1, o.data(), 1); };
    // returns norm of the vector
    double norm() { return std::sqrt(ddot(*this)); };
    // daxpy added to self
    void daxpy(double a, RotFile& o) { daxpy_(size_, a, o.data(), 1, data(), 1); }; 
    // orthogonalize to the liset of RotFile's
    double orthog(std::list<std::shared_ptr<RotFile> > c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter)
        this->daxpy(- this->ddot(**iter), **iter);
      const double scal = 1.0/this->norm();
      dscal_(size_, scal, data(), 1);
      return 1.0/scal;
    };

    // return data_
    double* data() { return data_.get(); };
    const double* data() const { return data_.get(); };
    // return data_
    double* begin() { return data(); };
    // return data_
    double* end() { return data()+size_; };

    // closed-active block. closed runs first
    double* ptr_ca() { return data(); };
    double& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; };
    // active-virtual block. virtual runs first
    double* ptr_va() { return data() + nclosed_*nact_; };
    double& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; };
    // closed-virtual block. virtual runs first
    double* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; };
    double& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; };
    // reference config.
    double& ele_ref() { assert(superci_); return data_[size_-1]; };
    // const references
    const double& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; };
    const double& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; };
    const double& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; };
    const double& ele_ref() const { assert(superci_); return data_[size_-1]; };

    // unpack to Matrix1e
    std::shared_ptr<Matrix1e> unpack(std::shared_ptr<Geometry> geom) const;

    // print matrix
    void print() const;
};


class QFile {
  protected:
    std::unique_ptr<double[]> data_;
    int na_;
    int nb_;

  public:
    QFile(const int nb, const int na) : data_(new double[na*nb]), na_(na), nb_(nb) {
      std::fill(data(), data()+na*nb, 0.0);
    };
    QFile(const QFile& o) : data_(new double[o.na_*o.nb_]), na_(o.na_), nb_(o.nb_) {
      std::copy(o.data(), o.data()+na_*nb_, data());
    };
    QFile(const std::shared_ptr<QFile> o) : data_(new double[o->na_*o->nb_]), na_(o->na_), nb_(o->nb_) {
      std::copy(o->data(), o->data()+na_*nb_, data());
    };
    ~QFile() { };

    double* data() { return data_.get(); };
    const double* data() const { return data_.get(); };
    double& element(const int i, const int j) { return data_[i+j*nb_]; };
    double* element_ptr(const int i, const int j) { return data()+i+j*nb_; };

    void zero() { std::fill(data(), data()+na_*nb_, 0.0); };

    QFile operator*(const QFile& o) {
      const int n = nb_;
      const int nc = na_; assert(nc == o.nb_);
      const int m = o.na_;
      QFile out(m, n);
      dgemm_("N", "N", n, m, nc, 1.0, data(), n, o.data(), nc, 0.0, out.data(), n); 
      return out;
    };

    QFile& operator*=(const double a) { dscal_(nb_*na_, a, data(), 1); };

    QFile operator*(const double a) {
      QFile tmp(*this);
      dscal_(nb_*na_, a, tmp.data(), 1);
      return tmp;
    };

    QFile operator/(const double a) {
      QFile tmp(*this);
      dscal_(nb_*na_, 1.0/a, tmp.data(), 1);
      return tmp;
    };

    QFile operator%(const QFile& o) {
      const int n = na_;
      const int nc = nb_; assert(nc == o.nb_);
      const int m = o.na_;
      QFile out(m, n);
      dgemm_("T", "N", n, m, nc, 1.0, data(), n, o.data(), nc, 0.0, out.data(), n); 
      return out;
    };

    QFile operator^(const QFile& o) {
      const int n = nb_;
      const int nc = na_; assert(nc == o.na_);
      const int m = o.nb_;
      QFile out(m, n);
      dgemm_("N", "T", n, m, nc, 1.0, data(), n, o.data(), nc, 0.0, out.data(), n); 
      return out;
    };


    QFile& operator+=(const QFile& o) {
      daxpy_(nb_*na_, 1.0, o.data(), 1, data(), 1); 
    };

    QFile operator+(const QFile& o) {
      QFile tmp(*this);
      daxpy_(nb_*na_, 1.0, o.data(), 1, tmp.data(), 1); 
      return tmp;
    };

    QFile operator-(const QFile& o) {
      QFile tmp(*this);
      daxpy_(nb_*na_, -1.0, o.data(), 1, tmp.data(), 1); 
      return tmp;
    };


    QFile& operator=(const QFile& o) {
      dcopy_(nb_*na_, o.data(), 1, data(), 1);
    };

    void print() const {
      std::cout << " -- printing QFile content --" << std::endl;
      const double* d = data();
      for (int i = 0; i != na_; ++i) {
        for (int j = 0; j != nb_; ++j, ++d)
          std::cout << std::setw(12) << std::setprecision(8) << *d; 
        std::cout << std::fixed << std::endl;
      }
    };
};

#endif
