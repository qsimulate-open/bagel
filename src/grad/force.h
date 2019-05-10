//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: force.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_GRAD_FORCE_H
#define __SRC_GRAD_FORCE_H

#include <src/wfn/reference.h>
#include <src/util/muffle.h>
#include <src/grad/gradinfo.h>

namespace bagel {

class Force {
  protected:
    const std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;

    bool numerical_;
    std::vector<double> energy_;
    std::vector<double> force_dipole_;

  public:
    Force(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    std::shared_ptr<GradFile> compute();
    void force_export(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo, const std::vector<double> energy, std::shared_ptr<const GradFile> out, const bool export_single);
    const std::vector<double>& force_dipole() const { return force_dipole_; }

    std::shared_ptr<const Reference> conv_to_ref() const { return ref_; }
};

}

#endif
