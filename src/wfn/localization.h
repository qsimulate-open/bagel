//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: localization.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __BAGEL_WFN_LOCALIZE_H
#define __BAGEL_WFN_LOCALIZE_H

#include <vector>

#include <src/wfn/reference.h>

namespace bagel {

class OrbitalLocalization {
  protected:
    // Hang on to the input for convenience
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Matrix> coeff_;

    std::vector<std::pair<int, int>> orbital_subspaces_;
    std::vector<std::pair<int, int>> region_bounds_;

    // used to reorder subspaces only
    VectorB diagonals_;

    virtual std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> coeff) = 0;

  public:
    OrbitalLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<std::pair<int, int>> subspaces);
    OrbitalLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref);

    std::vector<std::pair<int, int>> orbital_subspaces() const { return orbital_subspaces_; }

    std::shared_ptr<Matrix> localize();
    virtual double metric() const = 0;
};

class RegionLocalization : public OrbitalLocalization {
  protected:
    std::vector<int> sizes_;
    std::shared_ptr<Matrix> sqrt_S_;
    std::shared_ptr<Matrix> S_inverse_half_;

    std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> coeff) override;

  public:
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<std::pair<int, int>> subspaces, std::vector<int> region_sizes = std::vector<int>());
    RegionLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, std::vector<int> region_sizes = std::vector<int>());

    std::shared_ptr<Matrix> localize();

    double metric() const override { return 0.0; } // Until we can think of a good metric

  private:
    void common_init(std::vector<int> sizes);
};

// Pipek-Mezey
class PMLocalization : public OrbitalLocalization {
  protected:
    std::shared_ptr<Matrix> S_;

    int max_iter_;
    double thresh_;
    bool lowdin_;

    std::shared_ptr<Matrix> localize_space(std::shared_ptr<const Matrix> coeff) override;

  public:
    PMLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Matrix> coeff,
      std::vector<std::pair<int, int>> subspaces, std::vector<int> region_sizes = std::vector<int>());
    PMLocalization(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, std::vector<int> region_sizes = std::vector<int>());

    double metric() const override;

  private:
    double calc_P(std::shared_ptr<const Matrix> coeff, const int nstart, const int norb) const;
    void common_init(std::vector<int> sizes);
};

}

#endif
