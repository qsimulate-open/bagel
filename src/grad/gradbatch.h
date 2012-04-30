//
// Newint - Parallel electron correlation program.
// Filename: gradbatch.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#ifndef __SRC_RYSINT_GRAD_H
#define __SRC_RYSINT_GRAD_H

// compute analytic nuclear gradients
#include <vector>
#include <memory>
#include <src/scf/shell.h>
#include <src/rysint/eribatch.h>

class GradBatch : public RysInt {
  protected:
    // if we only compute three-center integrals, we want to use this info
    // to reduce the number of differentiation
    int centers_;

  public:
    GradBatch(const std::vector<std::shared_ptr<Shell> > shells, const double max_density, const double dummy = 0.0, const bool dum = true)
      : RysInt(shells) {
      centers_ = 4;  
      for (auto i = shells.begin(); i != shells.end(); ++i)
        if ((*i)->dummy()) --centers_;

      // a member variable in RysInt <- ERIBatch <- GradBatch.
      deriv_rank_ = 1;
    };
    ~GradBatch() {};

    void compute();

};

#endif

