//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexdf_base.h
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

// This base class contains utility functions needed by the various ComplexDF classes.
// It has only the minimal functionality not also needed for standard DF classes.


#ifndef __SRC_DF_COMPLEXDF_BASE_H
#define __SRC_DF_COMPLEXDF_BASE_H

#include <src/df/dfblock.h>
#include <src/df/paralleldf.h>

namespace bagel {

class ComplexDF_base {

  protected:
    // each one is assigned to half the members of ParallelDF::block_
    std::vector<std::shared_ptr<DFBlock>> real_block_;
    std::vector<std::shared_ptr<DFBlock>> imag_block_;

    // these would be better obtained from ParallelDF directly...
    int cnaux_;
    bool cserial_;

  public:
    ComplexDF_base() { }

    std::shared_ptr<const DFBlock> real_block(int i) const { return real_block_[i]; }
    std::shared_ptr<const DFBlock> imag_block(int i) const { return imag_block_[i]; }

    std::shared_ptr<ZVectorB> complex_compute_cd(const std::shared_ptr<const ZMatrix> den, std::shared_ptr<const Matrix> dat2, const bool onlyonce = false) const;

    // call in constructor of derived class
    void assign_complex_blocks(ParallelDF& source);

};

}

#endif
