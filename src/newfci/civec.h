//
// BAGEL - Parallel electron correlation program.
// Filename: civec.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef NEWINT_NEWFCI_CIVEC_H
#define NEWINT_NEWFCI_CIVEC_H

#include <cstddef>
#include <list>
#include <memory>
#include <vector>
#include <algorithm>
#include <cassert>
#include <src/util/f77.h>
#include <src/newfci/determinants.h>

namespace bagel {

class NewCivec {
  protected:
    // The determinant space in which this NewCivec object is defined
    mutable std::shared_ptr<const NewDeterminants> det_;

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
    const double* cc() const { return cc_ptr_; };

  public:
    NewCivec(std::shared_ptr<const NewDeterminants> det);

    // constructor that is called by Dvec.
    NewCivec(std::shared_ptr<const NewDeterminants> det, double* din_);

    // copy constructor
    NewCivec(const NewCivec& o);

    // this is not a copy constructor.
    NewCivec(std::shared_ptr<NewCivec> o, std::shared_ptr<const NewDeterminants> det);

    ~NewCivec() { };

    double* data() { return cc(); };
    double& element(size_t i, size_t j) { return cc(i+j*lenb_); }; // I RUNS FIRST
    double* element_ptr(size_t i, size_t j) { return cc()+i+j*lenb_; }; // I RUNS FIRST

    double& data(const int& i) { return cc(i); };
    const double& data(const int& i) const { return cc(i); };

    const double* data() const { return cc(); };
    const double* element_ptr(size_t i, size_t j) const { return cc()+i+j*lenb_; }; // I RUNS FIRST

    std::shared_ptr<const NewDeterminants> det() const { return det_; };
    void set_det(std::shared_ptr<const NewDeterminants> o) const { det_ = o; };

    void zero() { std::fill(cc(), cc()+lena_*lenb_, 0.0); };

    int size() const { return lena_*lenb_; };

    std::shared_ptr<NewCivec> transpose() const;

    int lena() const { return lena_; };
    int lenb() const { return lenb_; };

    // some functions for convenience
    void daxpy(double a, const NewCivec& other);
    double ddot(const NewCivec& other) const;
    double norm() const;
    double variance() const;
    void scale(const double a);

    NewCivec& operator*=(const double& a) { scale(a); return *this; };
    NewCivec& operator+=(const double& a) { daxpy_(size(),  1.0, &a, 0, data(), 1); return *this; }; // <- note I used a stride 0
    NewCivec& operator-=(const double& a) { daxpy_(size(), -1.0, &a, 0, data(), 1); return *this; };

    NewCivec& operator/=(const NewCivec& o);
    NewCivec operator/(const NewCivec& o) const;

    // assumes that NewCivec's in c are already orthogonal with each other.
    // returns scaling factor (see implementation)

    double orthog(std::list<std::shared_ptr<const NewCivec> > c);
    double orthog(std::shared_ptr<const NewCivec> o);
    void project_out(std::shared_ptr<const NewCivec> o) { daxpy(-ddot(*o), *o); }

    void print(const double& thresh) const { det_->print(data(), thresh); };
};

}

#endif
