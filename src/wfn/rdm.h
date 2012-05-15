//
// Newint - Parallel electron correlation program.
// Filename: rdm.h
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


#ifndef __NEWINT_WFN_RDM_H
#define __NEWINT_WFN_RDM_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <src/util/f77.h>

template <int rank>
class RDM {
  protected:
    std::unique_ptr<double[]> data_;
    const int norb_;
    int dim_;

  public:
    RDM(const int n) : norb_(n) {
      assert(rank > 0);
      dim_ = 1;
      for (int i = 0; i != rank; ++i) dim_ *= n;
      data_ = std::unique_ptr<double[]>(new double[dim_*dim_]);
    };
    RDM(const RDM& o) : norb_(o.norb_), dim_(o.dim_) {
      data_ = std::unique_ptr<double[]>(new double[dim_*dim_]);
      std::copy(o.data(), o.data()+dim_*dim_, data());
    };
    ~RDM() {  };

    double* data() { return data_.get(); };
    const double* data() const { return data_.get(); };
    double* first() { return data(); };
    double& element(int i, int j) { assert(i<dim_ && j<dim_); return data_[i+j*dim_]; };
    const double& element(int i, int j) const { assert(i<dim_ && j<dim_); return data_[i+j*dim_]; };
    const double* element_ptr(int i, int j) const { assert(i<dim_ && j<dim_); return &element(i,j); };
    // careful, this should not be called for those except for 2RDM.
    double& element(int i, int j, int k, int l) { assert(rank == 2 && i<norb_ && j<norb_ && k<norb_ && l<norb_);
        return data_[i+norb_*(j+norb_*(k+norb_*l))]; };
    double* element_ptr(int i, int j, int k, int l) { return &element(i,j,k,l); };
    const double& element(int i, int j, int k, int l) const { assert(rank == 2 && i<norb_ && j<norb_ && k<norb_ && l<norb_);
        return data_[i+norb_*(j+norb_*(k+norb_*l))]; };

    void zero() { std::fill(data(), data()+dim_*dim_, 0.0); };
    void daxpy(const double a, const RDM& o) {
      daxpy_(dim_*dim_, a, o.data(), 1, data(), 1);
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
      for (int i = 0; i != dim_; ++i) buf[i+i*dim_] = 2.0; 
      daxpy_(dim_*dim_, -1.0, data(), 1, &(buf[0]), 1);
      int lwork = 5*dim_;
      std::unique_ptr<double[]> work(new double[lwork]);
      int info;
      dsyev_("V", "U", &dim_, &(buf[0]), &dim_, &(vec[0]), work.get(), &lwork, &info);
      assert(!info);
      for (auto i = vec.begin(); i != vec.end(); ++i) *i = 2.0-*i;
      return std::make_pair(buf, vec);
    };

    void transform(const std::vector<double>& coeff) {
      const double* start = &(coeff[0]);
      std::unique_ptr<double[]> buf(new double[dim_*dim_]);
      if (rank == 1) {
        dgemm_("N", "N", dim_, dim_, dim_, 1.0, data(), dim_, start, dim_, 0.0, buf.get(), dim_); 
        dgemm_("T", "N", dim_, dim_, dim_, 1.0, start, dim_, buf.get(), dim_, 0.0, data(), dim_); 
      } else if (rank == 2) {
        // first half transformation
        dgemm_("N", "N", dim_*norb_, norb_, norb_, 1.0, data(), dim_*norb_, start, norb_, 0.0, buf.get(), dim_*norb_); 
        for (int i = 0; i != norb_; ++i)
          dgemm_("N", "N", dim_, norb_, norb_, 1.0, buf.get()+i*dim_*norb_, dim_, start, norb_, 0.0, data()+i*dim_*norb_, dim_); 
        // then tranpose
        mytranspose_(data(), &dim_, &dim_, buf.get());
        // and do it again
        dgemm_("N", "N", dim_*norb_, norb_, norb_, 1.0, buf.get(), dim_*norb_, start, norb_, 0.0, data(), dim_*norb_); 
        for (int i = 0; i != norb_; ++i)
          dgemm_("N", "N", dim_, norb_, norb_, 1.0, data()+i*dim_*norb_, dim_, start, norb_, 0.0, buf.get()+i*dim_*norb_, dim_);
        // to make sure for non-symmetric density matrices (and anyway this should be cheap).
        mytranspose_(buf.get(), &dim_, &dim_, data());
      } else {
        assert(false);
      }
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
