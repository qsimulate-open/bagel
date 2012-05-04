//
// Newint - Parallel electron correlation program.
// Filename: gnaibatch.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#ifndef __SRC_GRAD_GNAIBATCH_H
#define __SRC_GRAD_GNAIBATCH_H

#include <memory>
#include <src/rysint/naibatch_base.h>

class GNAIBatch : public NAIBatch_base {

  protected:
    void set_exponents() {
      exponents_ = std::unique_ptr<double[]>(new double[primsize_*2]);
      double* tmp = exponents_.get();
      for (auto i0 = basisinfo_[0]->exponents().begin(); i0 != basisinfo_[0]->exponents().end(); ++i0) {
        for (auto i1 = basisinfo_[1]->exponents().begin(); i1 != basisinfo_[1]->exponents().end(); ++i1, tmp+=2) {
          tmp[0] = *i0;
          tmp[1] = *i1;
        }
      }
    };
    std::unique_ptr<double[]> exponents_;

  public:
    
    GNAIBatch(const std::vector<std::shared_ptr<Shell> > _info, const std::shared_ptr<Geometry> gm, const int L = 0, const double A = 0.0)
      :  NAIBatch_base(_info, gm, 1, L, A) {
      set_exponents();
    };
    ~GNAIBatch() {};

    /// compute a batch of integrals
    void compute();

};

#endif

