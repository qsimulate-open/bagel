
// BAGEL - Parallel electron correlation program.
// Filename: properties.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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



#ifndef __BAGEL_FCI_PROPERTIES_H
#define __BAGEL_FCI_PROPERTIES_H

#include <src/wfn/reference.h>
#include <src/scf/scf.h>
#include <src/fci/dvec.h>

namespace bagel {

class CIProperties {
  protected:
    int nocc_;
    int norb_;
    int nbasis_;
    int sizeij_;

    double core_prop_;

    std::shared_ptr<const Matrix> prop_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const Coeff> coeff_;

  public:
    CIProperties(const std::shared_ptr<const Reference> ref, const int nstart, const int nfence)
      : nocc_(nstart), norb_(nfence-nstart), nbasis_(ref->geom()->nbasis()), geom_(ref->geom()), ref_(ref), coeff_(ref->coeff()) { }
    CIProperties(const std::shared_ptr<const Reference> ref, const int nstart, const int nfence, const std::shared_ptr<const Coeff> coeff)
      : nocc_(nstart), norb_(nfence-nstart), nbasis_(ref->geom()->nbasis()), geom_(ref->geom()), ref_(ref), coeff_(coeff) { }

    virtual void compute(std::shared_ptr<const Dvec>) = 0;

    const std::shared_ptr<const Geometry> geom() const { return geom_; };
    const std::shared_ptr<const Matrix> prop() const { return prop_; };

    virtual double core() const { return core_prop_; };
    virtual void print() const = 0;

  private:
    virtual void init(const int nstart, const int nfence) = 0; // init should do everything it can that does not require the Civecs
};

// 1e properties of a CI wavefunction : This is kind of useless right now, but I envision it being less useless when there are more properties to calculate
// It just so happens that the only property I'm calculating at the moment is a special case (with 3 components)
class Prop1e : public CIProperties {
  public:
    Prop1e(const std::shared_ptr<const Reference> ref, const int nstart, const int nfence) : CIProperties(ref, nstart, nfence) {}
    Prop1e(const std::shared_ptr<const Reference> ref, const int nstart, const int nfence, const std::shared_ptr<const Coeff> coeff)
      : CIProperties(ref,nstart,nfence, coeff) {}
};

// Dipole operator of a CI wavefunction
class CIDipole : public Prop1e {
  protected:
    std::array<std::shared_ptr<const Matrix>, 3> dipole_mo_;
    std::array<std::shared_ptr<Matrix>, 3> dipole_matrices_;
    std::array<double,3> core_dipole_;
    std::array<std::unique_ptr<double[]>,3> compressed_dipoles_;

  public:
    CIDipole(const std::shared_ptr<const Reference> ref, const int nstart, const int nfence) : Prop1e(ref, nstart, nfence)
      { init(nstart,nfence); }
    CIDipole(const std::shared_ptr<const Reference> ref, const int nstart, const int nfence, const std::shared_ptr<const Coeff> coeff)
      : Prop1e(ref, nstart, nfence, coeff) { init(nstart, nfence); }

    void init(const int nstart, const int nfence) override;
    virtual void compute(std::shared_ptr<const Dvec> ccvec) override;

    std::array<double,3> core_dipole() const { return core_dipole_; }
    double core_dipole(const int i) const { return core_dipole_[i]; }
    std::array<std::shared_ptr<Matrix>, 3> dipoles() const { return dipole_matrices_; }
    std::shared_ptr<Matrix> dipoles(const int i) const { return dipole_matrices_[i]; }

    virtual void print() const override;
};

}

#endif
