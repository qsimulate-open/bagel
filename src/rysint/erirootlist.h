//
// Newint - Parallel electron correlation program.
// Filename: erirootlist.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __src_rysint_erirootlist_h
#define __src_rysint_erirootlist_h

#include <src/rysint/macros.h>

struct ERIRootList  {
  public:
    ERIRootList();
    ~ERIRootList();

    void (*rfunc[RYS_MAX + 1])(const double*, double*, double*, const int*);

    void rootfunc_call(const unsigned int i, const double* a1, double* a2, double* a3, const int* a4) {
      return (rfunc[i])(a1, a2, a3, a4);
    };

}; 

#endif

