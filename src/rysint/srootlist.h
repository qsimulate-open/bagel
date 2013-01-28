//
// BAGEL - Parallel electron correlation program.
// Filename: srootlist.h
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


#ifndef __SRC_RYSINT_SROOTLIST_H
#define __SRC_RYSINT_SROOTLIST_H

#include <functional>
#include <src/rysint/intparam.h>

struct SRootList  {
  public:
    SRootList();
    ~SRootList();

    std::function<void (const double*, const double*, double*, double*, const int*)> srfunc[RYS_MAX + 1];

    void srootfunc_call(const unsigned int i, const double* a0, const double* a1, double* a2, double* a3, const int* a4) {
      return (srfunc[i])(a0, a1, a2, a3, a4);
    };

};

#endif

