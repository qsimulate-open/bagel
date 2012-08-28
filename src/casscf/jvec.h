//
// BAGEL - Parallel electron correlation program.
// Filename: jvec.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __NEWINT_SRC_CASSCF_JVEC_H
#define __NEWINT_SRC_CASSCF_JVEC_H

#include <src/fci/fci.h>
#include <src/scf/coeff.h>
#include <src/df/df.h>

namespace bagel {

class Jvec {
  protected:
    std::shared_ptr<const DF_Half> half_; 
    std::shared_ptr<const DF_Full> jvec_;
    std::unique_ptr<double[]> rdm2_all_; 

  public:
    Jvec(std::shared_ptr<FCI> fci, std::shared_ptr<const Coeff> c, const size_t, const size_t, const size_t);
    ~Jvec() {};

    const std::shared_ptr<const DF_Half> half() const { return half_; };
    const std::shared_ptr<const DF_Full> jvec() const { return jvec_; };
    const double* half_ptr() const { return half_->data(); };
    const double* jvec_ptr() const { return jvec_->data(); };
    const double* rdm2_all() const { return rdm2_all_.get(); };

};

}

#endif
