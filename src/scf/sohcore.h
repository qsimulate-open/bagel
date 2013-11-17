//
// BAGEL - Parallel electron correlation program.
// Filename: sohcore.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_SCF_SOHCORE_H
#define __SRC_SCF_SOHCORE_H

#include <src/math/matrix.h>
#include <src/wfn/geometry.h>
#include <src/scf/sohcore_base.h>

namespace bagel {

class SOHcore : public Matrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const SOHcore_base> hcore_;
//  std::shared_ptr<const Matrix> so1_;
//  std::shared_ptr<const Matrix> so2_;
//  std::shared_ptr<const Matrix> ecp_;
    std::shared_ptr<const Matrix> so1();
    std::shared_ptr<const Matrix> so2();
    std::shared_ptr<const Matrix> ecp();

    // initialize "this"
    void form_sohcore();

  public:
    SOHcore(const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const SOHcore_base> h);

};

}

#endif

