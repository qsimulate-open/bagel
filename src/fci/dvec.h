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


#ifndef BAGEL_FCI_DVEC_H
#define BAGEL_FCI_DVEC_H

#include <src/fci/civec.h>

// I forgot why I named this class "Dvec"...
// Basically Dvec is a vector of Civec, with some utility functions.
// It is used as intermediates in FCI, or eigenvectors in multistate CASSCF.

// I decided to make this class, in order to ensure that the date area
// in a FCI intermediate be consecutive so that I can use DGEMM etc.


// TODO The Dvec class is NOT yet flexible for Civectors with different Determinants objects.
// This can be easily done by modifing what follows.

namespace bagel {

class Dvec {
  protected:
    // the determinant space where Dvec's are sitting
    mutable std::shared_ptr<const Determinants> det_;

    size_t lena_;
    size_t lenb_;
    // the size of the vector<shared_ptr<Civec>>
    size_t ij_;
    std::vector<std::shared_ptr<Civec>> dvec_;
    std::unique_ptr<double[]> data_;

  public:
    Dvec(std::shared_ptr<const Determinants> det, const size_t ij);

    Dvec(const Dvec& o);

    Dvec(std::shared_ptr<const Civec> e, const size_t ij);

    // I think this is very confusiong... this is done this way in order not to delete Civec when Dvec is deleted.
    Dvec(std::shared_ptr<const Dvec> o);

    Dvec(std::vector<std::shared_ptr<Civec>> o);

    ~Dvec() { }

    std::shared_ptr<const Determinants> det() const { return det_; }

    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }

    std::shared_ptr<Civec>& data(const size_t i) { return dvec_[i]; }
    std::shared_ptr<const Civec> data(const size_t i) const { return dvec_[i]; }
    void zero() { std::fill(data(), data()+lena_*lenb_*ij_, 0.0); }

    std::vector<std::shared_ptr<Civec>>& dvec() { return dvec_; }
    const std::vector<std::shared_ptr<Civec>>& dvec() const { return dvec_; }
    std::vector<std::shared_ptr<Civec>> dvec(const std::vector<int>& conv);
    std::vector<std::shared_ptr<const Civec>> dvec(const std::vector<int>& conv) const;

    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }
    size_t ij() const { return ij_; }
    size_t size() const { return lena_*lenb_*ij_; }

    void set_det(std::shared_ptr<const Determinants> o) const;

    // some functions for convenience
    double ddot(const Dvec& other) const;
    void daxpy(double a, const Dvec& other);
    Dvec& operator+=(const Dvec& o) { daxpy(1.0, o); return *this; }
    Dvec& operator-=(const Dvec& o) { daxpy(-1.0, o); return *this; }

    Dvec operator+(const Dvec& o) const { Dvec out(*this); return out += o; }
    Dvec operator-(const Dvec& o) const { Dvec out(*this); return out -= o; }

    Dvec& operator/=(const Dvec& o);
    Dvec operator/(const Dvec& o) const;

    double norm() const { return std::sqrt(ddot(*this)); }
    void scale(const double a);
    Dvec& operator*=(const double& a) { scale(a); return *this; }

    std::shared_ptr<Dvec> clone() const { return std::make_shared<Dvec>(det_, ij_); }
    std::shared_ptr<Dvec> copy() const { return std::make_shared<Dvec>(*this); }
    std::shared_ptr<Dvec> spin() const;
    std::shared_ptr<Dvec> spinflip(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;
    std::shared_ptr<Dvec> spin_lower(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;
    std::shared_ptr<Dvec> spin_raise(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;

    void orthog(std::shared_ptr<const Dvec> o);
    void project_out(std::shared_ptr<const Dvec> o);

    void print(const double thresh = 0.05) const;
};

}

#endif
