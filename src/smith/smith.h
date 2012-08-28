//
// BAGEL - Parallel electron correlation program.
// Filename: smith.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


// compiles some input data for the smith routines.
#ifndef __SMIT_SMITH_H
#define __SMIT_SMITH_H

namespace bagel {
namespace SMITH {

class SMITH_info {
  protected:
    int maxiter_;
    double thresh_residual_;

  public:
    SMITH_info();
    ~SMITH_info();

    int maxiter() const { return maxiter_; };
    double thresh_residual() const { return thresh_residual_; };
};

}
}

#endif
