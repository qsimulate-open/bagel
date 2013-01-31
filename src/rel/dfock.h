//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.h
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


#ifndef __SRC_REL_DFOCK_H
#define __SRC_REL_DFOCK_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/wfn/geometry.h>
#include <src/util/zmatrix.h>
#include <src/df/df.h>
#include <src/rel/dfhalfcomplex.h>
#include <src/rel/relhcore.h>

namespace bagel {

class DFock : public ZMatrix {
  protected:
    std::array<std::shared_ptr<DFHalfComplex>, 2> large_half_;
    std::array<std::shared_ptr<DFHalfComplex>, 18> small_half_;
    //std::vector<std::shared_ptr<Matrix>> rsmall_data_;
    //std::vector<std::shared_ptr<Matrix>> ismall_data_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const RelHcore> hcore_;
    void two_electron_part(const std::array<std::shared_ptr<const ZMatrix>, 4> ocoeff, const bool rhf, const double scale_ex);
    void compute_half_complex(std::array<std::shared_ptr<const Matrix>, 4>, std::array<std::shared_ptr<const Matrix>, 4>, 
                              std::shared_ptr<const DFDist>, std::vector<std::shared_ptr<DFDist> >);

    std::array<std::shared_ptr<const ZMatrix>, 4> ocoeff_;

  public:
    DFock(const std::shared_ptr<const Geometry> a, 
          const std::shared_ptr<const RelHcore> h,
          const std::shared_ptr<const ZMatrix> coeff, const bool rhf = true, const double scale_ex = 1.0)
     : ZMatrix(*h), geom_(a), hcore_(h) {
       
       assert(geom_->nbasis()*4 == coeff->ndim());
       for (int i = 0; i != 4; ++i)
         ocoeff_[i] = coeff->get_submatrix(i*geom_->nbasis(), 0, geom_->nbasis(), coeff->mdim()); 
       two_electron_part(ocoeff_, rhf, scale_ex);
    };

//    std::shared_ptr<Reference> conv_to_ref() const override;

};

}

#endif
