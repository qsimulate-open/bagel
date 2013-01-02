//
// BAGEL - Parallel electron correlation program.
// Filename: breitrootlist.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_RYSINT_BREITROOTLIST_H
#define __SRC_RYSINT_BREITROOTLIST_H

#include <src/rysint/intparam.h>
#include <src/rysint/f77.h>

namespace bagel {

struct BreitRootList  {
  private:
    void (*rfunc[RYS_MAX + 1])(const double*, double*, double*, const int*);

  public:
    BreitRootList() {
      rfunc[1] = &breitroot1_;
      rfunc[2] = &breitroot2_;
      rfunc[3] = &breitroot3_;
      rfunc[4] = &breitroot4_;
      rfunc[5] = &breitroot5_;
      rfunc[6] = &breitroot6_;
      rfunc[7] = &breitroot7_;
      rfunc[8] = &breitroot8_;
      rfunc[9] = &breitroot9_;
      rfunc[10] = &breitroot10_;
      rfunc[11] = &breitroot11_;
      rfunc[12] = &breitroot12_;
      rfunc[13] = &breitroot13_;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { (rfunc[i])(a1, a2, a3, &a4); }

};

}

#endif

