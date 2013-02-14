//
// BAGEL - Parallel electron correlation program.
// Filename: spin2rootlist.h
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


#ifndef __SRC_RYSINT_SPIN2ROOTLIST_H
#define __SRC_RYSINT_SPIN2ROOTLIST_H

#include <functional>
#include <src/rysint/intparam.h>
#include <src/rysint/f77.h>

namespace bagel {

struct Spin2RootList  {
  private:
    std::function<void (const double*, double*, double*, const int*)> rfunc[RYS_MAX + 1];

  public:
    Spin2RootList() {
      rfunc[1] = &spin2root1_;
      rfunc[2] = &spin2root2_;
      rfunc[3] = &spin2root3_;
      rfunc[4] = &spin2root4_;
      rfunc[5] = &spin2root5_;
      rfunc[6] = &spin2root6_;
      rfunc[7] = &spin2root7_;
      rfunc[8] = &spin2root8_;
      rfunc[9] = &spin2root9_;
      rfunc[10] = &spin2root10_;
      rfunc[11] = &spin2root11_;
      rfunc[12] = &spin2root12_;
      rfunc[13] = &spin2root13_;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { rfunc[i](a1, a2, a3, &a4); }

};

}

#endif

