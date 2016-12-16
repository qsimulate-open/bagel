//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: current.h
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


#ifndef __SRC_PROP_CURRENT_H
#define __SRC_PROP_CURRENT_H

#include <src/wfn/method.h>

namespace bagel {

class Current : public Method {

  friend class CurrentTask;

  protected:
    bool relativistic_;
    bool paramagnetic_;
    bool diamagnetic_;
    size_t ngrid_;

    std::array<double,3> inc_size_;
    std::array<size_t,3> ngrid_dim_;

    std::shared_ptr<const ZMatrix> density_;

    std::vector<double> coords_;

    // size = 3*(ngrid_+1); last 3 entries give the total integrated current
    std::vector<std::complex<double>> currents_;

    void computepoint(const size_t pos);
    void print() const;

  public:
    Current(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re);

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; };

};

}

#endif
