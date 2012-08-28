//
// BAGEL - Parallel electron correlation program.
// Filename: dvec.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef NEWINT_NEWFCI_DVEC_H
#define NEWINT_NEWFCI_DVEC_H

#include <src/newfci/civec.h>

// I forgot why I named this class "NewDvec"...
// Basically NewDvec is a vector of NewCivec, with some utility functions.
// It is used as intermediates in FCI, or eigenvectors in multistate CASSCF.

// I decided to make this class, in order to ensure that the date area
// in a FCI intermediate be consecutive so that I can use DGEMM etc.


// TODO The NewDvec class is NOT yet flexible for NewCivectors with different NewDeterminants objects.
// This can be easily done by modifing what follows.

class NewDvec {
  protected:
    // the determinant space where NewDvec's are sitting
    mutable std::shared_ptr<const NewDeterminants> det_;

    size_t lena_;
    size_t lenb_;
    // the size of the vector<shared_ptr<NewCivec> >
    size_t ij_;
    std::vector<std::shared_ptr<NewCivec> > dvec_;
    std::unique_ptr<double[]> data_;

  public:
    NewDvec(std::shared_ptr<const NewDeterminants> det, const size_t ij);

    NewDvec(const NewDvec& o);

    NewDvec(std::shared_ptr<const NewCivec> e, const size_t ij);

    // I think this is very confusiong... this is done this way in order not to delete NewCivec when NewDvec is deleted. 
    NewDvec(std::shared_ptr<const NewDvec> o);

    NewDvec(std::vector<std::shared_ptr<NewCivec> > o);

    ~NewDvec() { };

    std::shared_ptr<const NewDeterminants> det() const { return det_; };

    double* data() { return data_.get(); };
    const double* data() const { return data_.get(); };

    std::shared_ptr<NewCivec>& data(const size_t i) { return dvec_[i]; };
    std::shared_ptr<const NewCivec> data(const size_t i) const { return dvec_[i]; };
    void zero() { std::fill(data(), data()+lena_*lenb_*ij_, 0.0); };

    std::vector<std::shared_ptr<NewCivec> >& dvec() { return dvec_; };
    const std::vector<std::shared_ptr<NewCivec> >& dvec() const { return dvec_; };
    std::vector<std::shared_ptr<NewCivec> > dvec(const std::vector<int>& conv);
    std::vector<std::shared_ptr<const NewCivec> > dvec(const std::vector<int>& conv) const;

    size_t lena() const { return lena_; };
    size_t lenb() const { return lenb_; };
    size_t ij() const { return ij_; };
    size_t size() const { return lena_*lenb_*ij_; };

    void set_det(std::shared_ptr<const NewDeterminants> o) const;

    // some functions for convenience
    double ddot(const NewDvec& other) const;
    void daxpy(double a, const NewDvec& other);
    NewDvec& operator+=(const NewDvec& o) { daxpy(1.0, o); return *this; };
    NewDvec& operator-=(const NewDvec& o) { daxpy(-1.0, o); return *this; };

    NewDvec operator+(const NewDvec& o) const { NewDvec out(*this); return out += o; };
    NewDvec operator-(const NewDvec& o) const { NewDvec out(*this); return out -= o; };

    NewDvec& operator/=(const NewDvec& o);
    NewDvec operator/(const NewDvec& o) const;

    double norm() const { return std::sqrt(ddot(*this)); };
    void scale(const double a);
    NewDvec& operator*=(const double& a) { scale(a); return *this; };

    std::shared_ptr<NewDvec> clone() const { return std::shared_ptr<NewDvec>(new NewDvec(det_, ij_)); };
    std::shared_ptr<NewDvec> copy() const { return std::shared_ptr<NewDvec>(new NewDvec(*this)); };

    void orthog(std::shared_ptr<const NewDvec> o);
    void project_out(std::shared_ptr<const NewDvec> o);

    void print(const double thresh = 0.05) const;
};


#endif
