//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moprint.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_PROP_MOPRINT_H
#define __SRC_PROP_MOPRINT_H

#include <src/wfn/method.h>

namespace bagel {

class MOPrint : public Method {

  friend class MOPrintTask;

  protected:
    bool is_density_;
    bool relativistic_;
    bool paired_;
    bool cubefile_;
    size_t ngrid_;
    size_t norb_;

    std::array<double,3> inc_size_;
    std::array<size_t,3> ngrid_dim_;
    std::vector<int> orbitals_;

    std::vector<std::shared_ptr<const ZMatrix>> density_;

    std::vector<double> coords_;

    // size = ngrid_+1; last entry gives the total integrated charge
    std::vector<double> points_;

    void computepoint(const size_t pos);
    void computefull();
    void print() const;

  public:
    MOPrint(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re,
        const bool is_density = false, std::shared_ptr<const ZMatrix> total_density = std::make_shared<const ZMatrix>());

    void compute() override;
    bool is_density() const { return is_density_; };

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; };

};

}

#endif
