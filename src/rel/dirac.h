//
// BAGEL - Parallel electron correlation program.
// Filename: dirac.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_REL_DIRAC_H
#define __SRC_REL_DIRAC_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/wfn/geometry.h>
#include <src/scf/kinetic.h>
#include <src/util/matrix.h>
#include <src/rel/smallnai.h>
#include <src/rel/relhcore.h>
#include <src/rel/reloverlap.h>
#include <src/util/input.h>

namespace bagel {

class Dirac {
  protected:
    const std::shared_ptr<const Geometry> geom_;

    int max_iter_;
    int diis_start_;
    double thresh_scf_;
    double energy_;

    std::shared_ptr<const RelHcore> hcore_;
    std::shared_ptr<const RelOverlap> overlap_;
    std::shared_ptr<const RelOverlap> s12_;

  public:
    Dirac(const std::multimap<std::string, std::string>& idata_, const std::shared_ptr<const Geometry> geom,
          const std::shared_ptr<const Reference> re = std::shared_ptr<const Reference>());

    void compute();

    std::shared_ptr<Reference> conv_to_ref() const;

    void print_eig(const std::unique_ptr<double[]>&);

};

}

#endif
