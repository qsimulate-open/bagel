//
// BAGEL - Parallel electron correlation program.
// Filename: zdvec.h
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#ifndef BAGEL_ZFCI_ZDVEC_H
#define BAGEL_ZFCI_ZDVEC_H

#include <src/zfci/zcivec.h>

// I forgot why I named this class "Dvec"...
// Basically Dvec is a vector of Civec, with some utility functions.
// It is used as intermediates in FCI, or eigenvectors in multistate CASSCF.

// I decided to make this class, in order to ensure that the date area
// in a FCI intermediate be consecutive so that I can use DGEMM etc.


// TODO The Dvec class is NOT yet flexible for Civectors with different Determinants objects.
// This can be easily done by modifing what follows.

namespace bagel {

class ZDvec {
  protected:
    // the determinant space where Dvec's are sitting
    mutable std::shared_ptr<const Determinants> det_;

    size_t lena_;
    size_t lenb_;
    // the size of the vector<shared_ptr<Civec>>
    size_t ij_;
    std::vector<std::shared_ptr<ZCivec>> dvec_;
    std::unique_ptr<std::complex<double>[]> data_;

  public:
    ZDvec(std::shared_ptr<const Determinants> det, const size_t ij);

    ZDvec(const ZDvec& o);

    ZDvec(std::shared_ptr<const ZCivec> e, const size_t ij);

    // I think this is very confusiong... this is done this way in order not to delete Civec when Dvec is deleted.
    ZDvec(std::shared_ptr<const ZDvec> o);

    ZDvec(std::vector<std::shared_ptr<ZCivec>> o);

    ~ZDvec() { }

    std::shared_ptr<const Determinants> det() const { return det_; }

    std::complex<double>* data() { return data_.get(); }
    const std::complex<double>* data() const { return data_.get(); }

    std::shared_ptr<ZCivec>& data(const size_t i) { return dvec_[i]; }
    std::shared_ptr<const ZCivec> data(const size_t i) const { return dvec_[i]; }
    void zero() { std::fill(data(), data()+lena_*lenb_*ij_, 0.0); }

    std::vector<std::shared_ptr<ZCivec>>& dvec() { return dvec_; }
    const std::vector<std::shared_ptr<ZCivec>>& dvec() const { return dvec_; }
    std::vector<std::shared_ptr<ZCivec>> dvec(const std::vector<int>& conv);
    std::vector<std::shared_ptr<const ZCivec>> dvec(const std::vector<int>& conv) const;

    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }
    size_t ij() const { return ij_; }
    size_t size() const { return lena_*lenb_*ij_; }

    void set_det(std::shared_ptr<const Determinants> o) const;

    // some functions for convenience
    std::complex<double> zdotc(const ZDvec& other) const;
    void zaxpy(std::complex<double> a, const ZDvec& other);
    ZDvec& operator+=(const ZDvec& o) { zaxpy(1.0, o); return *this; }
    ZDvec& operator-=(const ZDvec& o) { zaxpy(-1.0, o); return *this; }

    ZDvec operator+(const ZDvec& o) const { ZDvec out(*this); return out += o; }
    ZDvec operator-(const ZDvec& o) const { ZDvec out(*this); return out -= o; }

    ZDvec& operator/=(const ZDvec& o);
    ZDvec operator/(const ZDvec& o) const;

    double norm() const { const std::complex<double> n = zdotc(*this); assert(fabs(n.imag()) < 1.0e-8); return std::sqrt(n.real()); }
    void scale(const std::complex<double> a);
    ZDvec& operator*=(const std::complex<double>& a) { scale(a); return *this; }

    std::shared_ptr<ZDvec> clone() const { return std::make_shared<ZDvec>(det_, ij_); }
    std::shared_ptr<ZDvec> copy() const { return std::make_shared<ZDvec>(*this); }
    std::shared_ptr<ZDvec> spin() const;
    std::shared_ptr<ZDvec> spinflip(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;
    std::shared_ptr<ZDvec> spin_lower(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;
    std::shared_ptr<ZDvec> spin_raise(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;

    void orthog(std::shared_ptr<const ZDvec> o);
    void project_out(std::shared_ptr<const ZDvec> o);

    void print(const double thresh = 0.05) const;
};

}

#endif
