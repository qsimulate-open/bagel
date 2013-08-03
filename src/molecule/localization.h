//
// BAGEL - Parallel electron correlation program.
// Filename: localization.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __BAGEL_MOLECULE_LOCALIZE_H
#define __BAGEL_MOLECULE_LOCALIZE_H

#include <vector>

#include <src/wfn/reference.h>

namespace bagel {

class OrbitalLocalization {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Matrix> coeff_;

    // These are set when the method is constructed
    const int nclosed_;
    const int nact_;
    const int nvirt_;

    // Can be toggled from the input
    bool localize_closed_;
    bool localize_active_;
    bool localize_virtual_;

    int max_iter_;
    double thresh_;

  public:
    OrbitalLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      const int nclosed, const int nact, const int nvirt);

    virtual std::shared_ptr<const Matrix> localize() = 0;
    virtual double metric() const = 0;
};

class RegionLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int>> bounds_;
    std::vector<int> sizes_;
    std::shared_ptr<Matrix> sqrt_S_;
    std::shared_ptr<Matrix> S_inverse_half_;
    std::vector<std::vector<int>> region_orbitals_;

  public:
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<int> region_sizes, const int nclosed, const int nact = 0, const int nvirt = 0);
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      const int nclosed, const int nact = 0, const int nvirt = 0);
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, std::vector<int> region_sizes) :
      RegionLocalization(input, ref->geom(), ref->coeff(), region_sizes, ref->nclosed(), ref->nact(), ref->nvirt()) {}
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref) :
      RegionLocalization(input, ref->geom(), ref->coeff(), ref->nclosed(), ref->nact(), ref->nvirt()) {}

    // TODO get rid of these or clean them up
    std::vector<std::vector<int>> region_orbitals() const { return region_orbitals_; }
    std::vector<int> region_orbitals(const int i) const { return region_orbitals_.at(i); }

    std::shared_ptr<const Matrix> localize() override;

    double metric() const override {return 0.0;} // Until we can think of a good metric

  private:
    void common_init(std::vector<int> sizes);
    std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> density);
};

// Pipek-Mezey
class PMLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int>> atom_bounds_;
    std::shared_ptr<Matrix> S_;

  public:
    PMLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      const int nclosed, const int nact = 0, const int nvirt = 0);
    PMLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref) :
      PMLocalization(input, ref->geom(), ref->coeff(), ref->nclosed(), ref->nact(), ref->nvirt()) {}

    std::shared_ptr<const Matrix> localize() override;

    double metric() const override;

  private:
    double calc_P(std::shared_ptr<const Matrix> coeff, const int nstart, const int norb) const;
    void common_init();

    void localize_space(std::shared_ptr<Matrix> coeff, const int nstart, const int norb);
};

}

#endif
