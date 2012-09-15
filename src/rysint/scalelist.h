//
// BAGEL - Parallel electron correlation program.
// Filename: scalelist.h
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


#ifndef __src_rysint_scalelist_h
#define __src_rysint_scalelist_h

#include <src/rysint/macros.h>

namespace bagel {

struct ScaleList {
  public:
    ScaleList();
    ~ScaleList();

    static void scale_data_1(double*, const double*, const double, const double*, const int);
    static void scale_data_2(double*, const double*, const double, const double*, const int);
    static void scale_data_3(double*, const double*, const double, const double*, const int);
    static void scale_data_4(double*, const double*, const double, const double*, const int);
    static void scale_data_5(double*, const double*, const double, const double*, const int);
    static void scale_data_6(double*, const double*, const double, const double*, const int);
    static void scale_data_7(double*, const double*, const double, const double*, const int);
    static void scale_data_8(double*, const double*, const double, const double*, const int);
    static void scale_data_9(double*, const double*, const double, const double*, const int);
    static void scale_data_10(double*, const double*, const double, const double*, const int);
    static void scale_data_11(double*, const double*, const double, const double*, const int);
    static void scale_data_12(double*, const double*, const double, const double*, const int);
    static void scale_data_13(double*, const double*, const double, const double*, const int);

    void (*scalefunc[RYS_MAX + 1])(double*, const double*, const double, const double*, const int);
};

}

#endif

