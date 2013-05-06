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


#ifndef __BAGEL_UTIL_LOCALIZE_H
#define __BAGEL_UTIL_LOCALIZE_H

#include <src/wfn/reference.h>

namespace bagel {

class OrbitalLocalization {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Coeff> coeff_;
    std::shared_ptr<const Reference> ref_;

    const int nclosed_;
    const int nact_;
    const int nvirt_;

    int iter_;
    double thresh_;

  public:
    OrbitalLocalization(std::shared_ptr<const Geometry> geom, std::shared_ptr<const Coeff> coeff, const int nclosed, const int nact = 0, const int nvirt = 0) : 
      geom_(geom), coeff_(coeff), nclosed_(nclosed), nact_(nact), nvirt_(nvirt) {}
    OrbitalLocalization(std::shared_ptr<const Reference> ref) : 
      OrbitalLocalization( ref->geom(), ref->coeff(), ref->nclosed(), ref->nact() ) { ref_ = ref; }

    virtual std::shared_ptr<const Coeff> localize(const int iter = 50, const double thresh = 1.0e-12) = 0;
    virtual double metric() const = 0;
};

class RegionLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int>> bounds_;
    std::vector<int> sizes_;
    std::shared_ptr<Matrix> sqrt_S_;
    std::shared_ptr<Matrix> S_inverse_half_;

  public:
    RegionLocalization(std::shared_ptr<const Geometry> geom, std::shared_ptr<const Coeff> coeff, std::vector<int> region_sizes, 
      const int nclosed, const int nact, const int nvirt = 0) : OrbitalLocalization(geom, coeff, nclosed, nact, nvirt) {common_init(region_sizes);}
    RegionLocalization(std::shared_ptr<const Reference> ref, std::vector<int> region_sizes) : 
      OrbitalLocalization(ref) {common_init(region_sizes);}

    std::shared_ptr<const Coeff> localize(const int iter = 0, const double thresh = 1.0e-12) override;

    double metric() const override {return 0.0;}
  
  private:
    void common_init(std::vector<int> sizes);
    std::shared_ptr<Matrix> localize_space(std::shared_ptr<Matrix> density);
};

// Pipek-Mezey
class PMLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int>> atom_bounds_;

    std::shared_ptr<Matrix> S_;

  public:
    PMLocalization(std::shared_ptr<const Geometry> geom, std::shared_ptr<const Coeff> coeff, const int nclosed, const int nact = 0, const int nvirt = 0) :
      OrbitalLocalization(geom, coeff, nclosed, nact, nvirt) {common_init(geom);}
    PMLocalization(std::shared_ptr<const Reference> ref) : OrbitalLocalization(ref) {common_init(ref->geom());}

    std::shared_ptr<const Coeff> localize(const int iter = 50, const double thresh = 1.0e-12) override;

    double metric() const override;

  private:
    double calc_P(std::shared_ptr<const Matrix> coeff, const int nstart, const int norb) const;
    void common_init(std::shared_ptr<const Geometry> geom);

    void localize_space(std::shared_ptr<Matrix> coeff, const int nstart, const int norb);
};

}

#endif
