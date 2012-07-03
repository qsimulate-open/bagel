//
// Newint - Parallel electron correlation program.
// Filename: civec.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef NEWINT_FCI_CIVEC_H
#define NEWINT_FCI_CIVEC_H

#include <list>
#include <memory>
#include <vector>
#include <algorithm>
#include <src/util/f77.h>

class Civec {
  protected:
    int lena_;
    int lenb_;

    // !!CAUTION!!
    // cc is formated so that B runs first.
    // Also, cc_ can be null if this is constructed by Dvec.
    std::unique_ptr<double[]> cc_;

    double* cc_ptr_;

    double& cc(int i) { return *(cc_ptr_+i); };
    const double& cc(int i) const { return *(cc_ptr_+i); };
    double* cc() { return cc_ptr_; };
    const double* const cc() const { return cc_ptr_; };

  public:
    Civec(const size_t lb, const size_t la) : lena_(la), lenb_(lb) {
      std::unique_ptr<double[]> cc__(new double[la*lb]);
      cc_ = std::move(cc__);
      cc_ptr_ = cc_.get();
      std::fill(cc(), cc() + la*lb, 0.0);
    }

    // constructor that is called by Dvec.
    Civec(const size_t lb, const size_t la, double* din_) : lena_(la), lenb_(lb) {
      cc_ptr_ = din_;
      std::fill(cc(), cc() + la*lb, 0.0);
    }

    // copy constructor
    Civec(const Civec& o) : lena_(o.lena()), lenb_(o.lenb()) {
      std::unique_ptr<double[]> cc__(new double[lena_*lenb_]);
      cc_ = std::move(cc__);
      cc_ptr_ = cc_.get();
      std::copy(o.cc(), o.cc() + lena_*lenb_, cc());
    };
    ~Civec() { };

    double& element(size_t i, size_t j) { return cc(i+j*lenb_); }; // I RUNS FIRST 
    double* element_ptr(size_t i, size_t j) { return cc()+i+j*lenb_; }; // I RUNS FIRST 
    double* data() { return cc(); };
    const double* data() const { return cc(); };


    void zero() { std::fill(cc(), cc()+lena_*lenb_, 0.0); };

    int size() const { return lena_*lenb_; };

    std::shared_ptr<Civec> transpose() {
      std::shared_ptr<Civec> ct(new Civec(lena_, lenb_));
      double* cct = ct->data(); 
      mytranspose_(cc(), &lenb_, &lena_, cct); 
      return ct;
    };
    int lena() const { return lena_; };
    int lenb() const { return lenb_; };

    // some functions for convenience
    double ddot(const Civec& other) const {
      return ddot_(lena_*lenb_, cc(), 1, other.data(), 1);
    };
    void daxpy(double a, const Civec& other) {
      daxpy_(lena_*lenb_, a, other.data(), 1, cc(), 1);
    };
    double norm() const {
      return std::sqrt(ddot_(lena_*lenb_, cc(), 1, cc(), 1));
    };
    void scale(const double a) {
      dscal_(lena_*lenb_, a, cc(), 1);
    };
    double variance() const {
      const int lab = lena_ * lenb_;
      return ddot_(lab, cc(), 1, cc(), 1)/lab;
    };

    // assumes that c is already orthogonal with each other.
    double orthog(std::list<std::shared_ptr<const Civec> > c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter) {
        const double scal = - this->ddot(**iter);
        this->daxpy(scal, **iter);
      }
      const double scal = 1.0/this->norm();
      dscal_(lena_*lenb_, scal, cc(), 1);
      return 1.0/scal; 
    };
};


class Dvec {
  protected:
    size_t lena_;
    size_t lenb_;
    size_t ij_;
    std::vector<std::shared_ptr<Civec> > dvec_;
    std::unique_ptr<double[]> data_;

  public:
    Dvec(const size_t lb, const size_t la, const size_t ij) : lena_(la), lenb_(lb), ij_(ij) {
      // actually data should be in a consecutive area to call dgemm.
      data_ = std::unique_ptr<double[]>(new double[lb*la*ij]);
      double* tmp = data();
      for (int i = 0; i != ij; ++i, tmp+=lb*la) {
        std::shared_ptr<Civec> c(new Civec(lb, la, tmp)); 
        dvec_.push_back(c);
      }
    };
    Dvec(const Dvec& o) : lena_(o.lena_), lenb_(o.lenb_), ij_(o.ij_) {
      data_ = std::unique_ptr<double[]>(new double[lena_*lenb_*ij_]);
      double* tmp = data_.get();
      for (int i = 0; i != ij_; ++i, tmp+=lenb_*lena_) {
        std::shared_ptr<Civec> c(new Civec(lenb_, lena_, tmp)); 
        dvec_.push_back(c);
      }
      std::copy(o.data(), o.data()+lena_*lenb_*ij_, data());
    };
    Dvec(std::shared_ptr<const Civec> e, const size_t ij) : lena_(e->lena()), lenb_(e->lenb()), ij_(ij) {
      // actually data should be in a consecutive area to call dgemm.
      data_ = std::unique_ptr<double[]>(new double[lena_*lenb_*ij]);
      double* tmp = data();
      for (int i = 0; i != ij; ++i, tmp+=lenb_*lena_) {
        std::shared_ptr<Civec> c(new Civec(lenb_, lena_, tmp)); 
        dcopy_(lenb_*lena_, e->data(), 1, c->data(), 1);
        dvec_.push_back(c);
      }
    };

    // I think this is very confusiong... this is done this way in order not to delete Civec when Dvec is deleted. 
    Dvec(std::shared_ptr<const Dvec> o) : lena_(o->lena_), lenb_(o->lenb_), ij_(o->ij_) {
      for (int i = 0; i != ij_; ++i) {
        std::shared_ptr<Civec> c(new Civec(*(o->data(i))));
        dvec_.push_back(c);
      }
    };

    Dvec(std::vector<std::shared_ptr<Civec> > o) : lena_(o.front()->lena()), lenb_(o.front()->lenb()), ij_(o.size()) {
      dvec_ = o;
    };

    ~Dvec() { };

    double* data() { return data_.get(); };
    const double* const data() const { return data_.get(); };

    std::shared_ptr<Civec>& data(const size_t i) { return dvec_[i]; };
    std::shared_ptr<const Civec> data(const size_t i) const { return std::const_pointer_cast<Civec>(dvec_[i]); };
    void zero() { std::fill(data(), data()+lena_*lenb_*ij_, 0.0); };

    std::vector<std::shared_ptr<Civec> > dvec() { return dvec_; };
    std::vector<std::shared_ptr<Civec> > dvec(const std::vector<int>& conv) {
      std::vector<std::shared_ptr<Civec> > out;
      int i = 0;
      for (auto iter = dvec_.begin(); iter != dvec_.end(); ++iter, ++i) {
        if (conv[i] == 0) out.push_back(*iter);
      }
      return out;
    };
    std::vector<std::shared_ptr<const Civec> > dvec(const std::vector<int>& conv) const {
      std::vector<std::shared_ptr<const Civec> > out;
      int i = 0;
      for (auto iter = dvec_.begin(); iter != dvec_.end(); ++iter, ++i) {
        if (conv[i] == 0) out.push_back(*iter);
      }
      return out;
    };

    size_t lena() const { return lena_; };
    size_t lenb() const { return lenb_; };
    size_t ij() const { return ij_; };

    // some functions for convenience
    double ddot(const Dvec& other) const {
      assert(ij() == other.ij());
      double sum = 0.0;
      for (auto i = dvec_.begin(), j = other.dvec_.begin(); i != dvec_.end(); ++i, ++j)
        sum += ddot_(lena_*lenb_, (*i)->data(), 1, (*j)->data(), 1);
      return sum;
    };
    void daxpy(double a, const Dvec& other) {
      auto j = other.dvec_.begin();
      for (auto i = dvec_.begin(); i != dvec_.end(); ++i, ++j)
        daxpy_(lena_*lenb_, a, (*j)->data(), 1, (*i)->data(), 1);
    };
    double norm() const { return std::sqrt(ddot(*this)); };
    void scale(const double a) {
      for (auto i = dvec_.begin(); i != dvec_.end(); ++i)
        dscal_(lena_*lenb_, a, (*i)->data(), 1);
    };

    std::shared_ptr<Dvec> clone() const { return std::shared_ptr<Dvec>(new Dvec(lenb(), lena(), ij())); };
    std::shared_ptr<Dvec> copy() const { return std::shared_ptr<Dvec>(new Dvec(*this)); };
};


#endif
