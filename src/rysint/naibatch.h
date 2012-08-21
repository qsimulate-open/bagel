//
// Newint - Parallel electron correlation program.
// Filename: naibatch.h
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

#ifndef __SRC_RYSINT_NAIBATCH_H
#define __SRC_RYSINT_NAIBATCH_H

#include <src/rysint/naibatch_base.h>

class NAIBatch : public NAIBatch_base {

  protected:

  public:
    
    NAIBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Geometry> gm, const int L = 0, const double A = 0.0)
      :  NAIBatch_base(_info, gm, 0, L, A) {};
     ~NAIBatch() {};

    /// compute a batch of integrals
    void compute();

};

#endif

