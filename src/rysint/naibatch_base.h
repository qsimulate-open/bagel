//
// Newint - Parallel electron correlation program.
// Filename: naibatch_base.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_RYSINT_NAIBATCH_BASE_H
#define __SRC_RYSINT_NAIBATCH_BASE_H

#include <memory>
#include <vector>
#include <tuple>
#include <src/rysint/macros.h>
#include <src/scf/shell.h>
#include <src/rysint/rysint.h>
#include <src/scf/geometry.h>


class NAIBatch_base : public RysInt {
  protected:
    std::shared_ptr<const Geometry> geom_;
    int natom_;

    /// for periodic calculations (UNCHECKED!!)
    const int L_;
    const double A_;

    void root_weight(const int ps); 
    void compute_ssss(const double);

  public:
    NAIBatch_base(const std::vector<std::shared_ptr<Shell> >& _info, const std::shared_ptr<const Geometry> gm, const int deriv,
                  const int L = 0, const double A = 0.0);
    ~NAIBatch_base() {};

    const std::shared_ptr<const Geometry> geom() const { return geom_; };


}; 

#endif
