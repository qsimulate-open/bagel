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


#ifndef __BAGEL_MOLECULE_LOCALIZE_H
#define __BAGEL_MOLECULE_LOCALIZE_H

#include <vector>

#include <src/wfn/reference.h>

namespace bagel {

class OrbitalLocalization {
  protected:
    // Hang on to the input for convenience
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Matrix> coeff_;

    std::vector<std::pair<int, int>> subspaces_;

    // used to reorder subspaces only
    std::vector<double> diagonals_;

    virtual std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> coeff) = 0;

  public:
    OrbitalLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<std::pair<int, int>> subspaces);
    OrbitalLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref);

    std::shared_ptr<Matrix> localize();
    virtual double metric() const = 0;
};

class RegionLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int>> bounds_;
    std::vector<int> sizes_;
    std::shared_ptr<Matrix> sqrt_S_;
    std::shared_ptr<DistMatrix> S_inverse_half_;
    std::vector<std::vector<int>> region_orbitals_;

    std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> coeff) override;

  public:
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<std::pair<int, int>> subspaces, std::vector<int> region_sizes = std::vector<int>());
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, std::vector<int> region_sizes = std::vector<int>());

    std::shared_ptr<Matrix> localize();

    double metric() const override {return 0.0;} // Until we can think of a good metric

  private:
    void common_init(std::vector<int> sizes);
};

// Pipek-Mezey
class PMLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int>> atom_bounds_;
    std::shared_ptr<Matrix> S_;

    int max_iter_;
    double thresh_;

    std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> coeff) override;

  public:
    PMLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<std::pair<int, int>> subspaces);
    PMLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref);

    double metric() const override;

  private:
    double calc_P(std::shared_ptr<const Matrix> coeff, const int nstart, const int norb) const;
    void common_init();
};

}

#endif
