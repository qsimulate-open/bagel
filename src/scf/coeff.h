//
// BAGEL - Parallel electron correlation program.
// Filename: coeff.h
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


#ifndef __src_scf_coeff_h
#define __src_scf_coeff_h

#include <src/scf/geometry.h>
#include <src/util/matrix.h>
#include <memory>

namespace bagel {

class Coeff : public Matrix {
  protected:
    std::shared_ptr<const Geometry> geom_;

  private:
//  std::shared_ptr<const Geometry> supergeom(std::vector<std::shared_ptr<const Coeff>> coeff_vec);
    int num_basis(std::vector<std::shared_ptr<const Coeff>> coeff_vec) const;

  public:
    Coeff(const Matrix&);
    Coeff(std::vector<std::shared_ptr<const Coeff>> coeff_vec);
    Coeff(std::shared_ptr<const Geometry> g) : Matrix(g->nbasis(), g->nbasis()), geom_(g) {};
    ~Coeff();

    std::shared_ptr<const Geometry> geom() const { assert(geom_); return geom_; };

    std::shared_ptr<Matrix> form_density_rhf(const int n, const int offset = 0) const;
    std::shared_ptr<Matrix> form_weighted_density_rhf(const int n, const std::vector<double>& e, const int offset = 0) const;
    std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> split(const int, const int) const;
};

}

#endif
