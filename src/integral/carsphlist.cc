//
// BAGEL - Parallel electron correlation program.
// Filename: carsphlist.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/integral/carsphlist.h>

// same convention as HRR

using namespace bagel;

CarSphList::CarSphList() {
      carsphfunc[ANG_HRR_END * 0 + 0] = &carsph_00;
      carsphfunc[ANG_HRR_END * 1 + 0] = &carsph_10;
      carsphfunc[ANG_HRR_END * 1 + 1] = &carsph_11;
      carsphfunc[ANG_HRR_END * 2 + 0] = &carsph_20;
      carsphfunc[ANG_HRR_END * 2 + 1] = &carsph_21;
      carsphfunc[ANG_HRR_END * 2 + 2] = &carsph_22;
      carsphfunc[ANG_HRR_END * 3 + 0] = &carsph_30;
      carsphfunc[ANG_HRR_END * 3 + 1] = &carsph_31;
      carsphfunc[ANG_HRR_END * 3 + 2] = &carsph_32;
      carsphfunc[ANG_HRR_END * 3 + 3] = &carsph_33;
      carsphfunc[ANG_HRR_END * 4 + 0] = &carsph_40;
      carsphfunc[ANG_HRR_END * 4 + 1] = &carsph_41;
      carsphfunc[ANG_HRR_END * 4 + 2] = &carsph_42;
      carsphfunc[ANG_HRR_END * 4 + 3] = &carsph_43;
      carsphfunc[ANG_HRR_END * 4 + 4] = &carsph_44;
      carsphfunc[ANG_HRR_END * 5 + 0] = &carsph_50;
      carsphfunc[ANG_HRR_END * 5 + 1] = &carsph_51;
      carsphfunc[ANG_HRR_END * 5 + 2] = &carsph_52;
      carsphfunc[ANG_HRR_END * 5 + 3] = &carsph_53;
      carsphfunc[ANG_HRR_END * 5 + 4] = &carsph_54;
      carsphfunc[ANG_HRR_END * 5 + 5] = &carsph_55;
      carsphfunc[ANG_HRR_END * 6 + 0] = &carsph_60;
      carsphfunc[ANG_HRR_END * 6 + 1] = &carsph_61;
      carsphfunc[ANG_HRR_END * 6 + 2] = &carsph_62;
      carsphfunc[ANG_HRR_END * 6 + 3] = &carsph_63;
      carsphfunc[ANG_HRR_END * 6 + 4] = &carsph_64;
      carsphfunc[ANG_HRR_END * 6 + 5] = &carsph_65;
      carsphfunc[ANG_HRR_END * 6 + 6] = &carsph_66;
}


CCarSphList::CCarSphList() {
      carsphfunc[ANG_HRR_END * 0 + 0] = &carsph_00;
      carsphfunc[ANG_HRR_END * 1 + 0] = &carsph_10;
      carsphfunc[ANG_HRR_END * 1 + 1] = &carsph_11;
      carsphfunc[ANG_HRR_END * 2 + 0] = &carsph_20;
      carsphfunc[ANG_HRR_END * 2 + 1] = &carsph_21;
      carsphfunc[ANG_HRR_END * 2 + 2] = &carsph_22;
      carsphfunc[ANG_HRR_END * 3 + 0] = &carsph_30;
      carsphfunc[ANG_HRR_END * 3 + 1] = &carsph_31;
      carsphfunc[ANG_HRR_END * 3 + 2] = &carsph_32;
      carsphfunc[ANG_HRR_END * 3 + 3] = &carsph_33;
      carsphfunc[ANG_HRR_END * 4 + 0] = &carsph_40;
      carsphfunc[ANG_HRR_END * 4 + 1] = &carsph_41;
      carsphfunc[ANG_HRR_END * 4 + 2] = &carsph_42;
      carsphfunc[ANG_HRR_END * 4 + 3] = &carsph_43;
      carsphfunc[ANG_HRR_END * 4 + 4] = &carsph_44;
      carsphfunc[ANG_HRR_END * 5 + 0] = &carsph_50;
      carsphfunc[ANG_HRR_END * 5 + 1] = &carsph_51;
      carsphfunc[ANG_HRR_END * 5 + 2] = &carsph_52;
      carsphfunc[ANG_HRR_END * 5 + 3] = &carsph_53;
      carsphfunc[ANG_HRR_END * 5 + 4] = &carsph_54;
      carsphfunc[ANG_HRR_END * 5 + 5] = &carsph_55;
      carsphfunc[ANG_HRR_END * 6 + 0] = &carsph_60;
      carsphfunc[ANG_HRR_END * 6 + 1] = &carsph_61;
      carsphfunc[ANG_HRR_END * 6 + 2] = &carsph_62;
      carsphfunc[ANG_HRR_END * 6 + 3] = &carsph_63;
      carsphfunc[ANG_HRR_END * 6 + 4] = &carsph_64;
      carsphfunc[ANG_HRR_END * 6 + 5] = &carsph_65;
      carsphfunc[ANG_HRR_END * 6 + 6] = &carsph_66;
}
