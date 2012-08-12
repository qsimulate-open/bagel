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
#include <tuple>
#include <src/rysint/naibatch_base.h>

class GNAIBatch : public NAIBatch_base {

  protected:
    void set_exponents();
    std::unique_ptr<double[]> exponents_;

    std::tuple<int,int> iatom_;

  public:
    
    GNAIBatch(const std::vector<std::shared_ptr<const Shell> > _info, const std::shared_ptr<const Geometry> gm, const std::tuple<int,int> i,
              const int L = 0, const double A = 0.0)
      :  NAIBatch_base(_info, gm, 1, L, A), iatom_(i) {
      if (swap01_) {
        std::swap(std::get<0>(iatom_), std::get<1>(iatom_));
      }
      set_exponents();
    };
    ~GNAIBatch() {};

    /// compute a batch of integrals
    void compute();

};

#endif

