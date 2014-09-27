//
// BAGEL - Parallel electron correlation program.
// Filename: current.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_LONDON_CURRENT_H
#define __SRC_LONDON_CURRENT_H

#include <src/wfn/method.h>

namespace bagel {

class Current : public Method {

  friend class CurrentTask;

  protected:
    bool relativistic_;
    bool paramagnetic_;
    bool diamagnetic_;
    size_t ngrid_;

    std::shared_ptr<const ZMatrix> density_;

    std::vector<double> coords_;

    // size = 3*(ngrid_+1); last 3 entries give the total integrated current
    std::vector<double> currents_;

    void computepoint(const size_t pos);
    void print() const;

  public:
    Current(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re);

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; };

};

}

#endif
