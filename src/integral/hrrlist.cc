//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hrrlist.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/integral/hrrlist.h>

using namespace bagel;

HRRList::HRRList() {
      hrrfunc[ANG_HRR_END * 1 + 1] = &perform_HRR_20_11;
      hrrfunc[ANG_HRR_END * 2 + 1] = &perform_HRR_30_21;
      hrrfunc[ANG_HRR_END * 2 + 2] = &perform_HRR_40_22;
      hrrfunc[ANG_HRR_END * 3 + 1] = &perform_HRR_40_31;
      hrrfunc[ANG_HRR_END * 3 + 2] = &perform_HRR_50_32;
      hrrfunc[ANG_HRR_END * 3 + 3] = &perform_HRR_60_33;
      hrrfunc[ANG_HRR_END * 4 + 1] = &perform_HRR_50_41;
      hrrfunc[ANG_HRR_END * 4 + 2] = &perform_HRR_60_42;
      hrrfunc[ANG_HRR_END * 4 + 3] = &perform_HRR_70_43;
      hrrfunc[ANG_HRR_END * 4 + 4] = &perform_HRR_80_44;
      hrrfunc[ANG_HRR_END * 5 + 1] = &perform_HRR_60_51;
      hrrfunc[ANG_HRR_END * 5 + 2] = &perform_HRR_70_52;
      hrrfunc[ANG_HRR_END * 5 + 3] = &perform_HRR_80_53;
      hrrfunc[ANG_HRR_END * 5 + 4] = &perform_HRR_90_54;
      hrrfunc[ANG_HRR_END * 5 + 5] = &perform_HRR_a0_55;
      hrrfunc[ANG_HRR_END * 6 + 1] = &perform_HRR_70_61;
      hrrfunc[ANG_HRR_END * 6 + 2] = &perform_HRR_80_62;
      hrrfunc[ANG_HRR_END * 6 + 3] = &perform_HRR_90_63;
      hrrfunc[ANG_HRR_END * 6 + 4] = &perform_HRR_a0_64;
      hrrfunc[ANG_HRR_END * 6 + 5] = &perform_HRR_b0_65;
      hrrfunc[ANG_HRR_END * 6 + 6] = &perform_HRR_c0_66;
#ifdef COMPILE_J_ORB
      hrrfunc[ANG_HRR_END * 7 + 1] = &perform_HRR_80_71;
      hrrfunc[ANG_HRR_END * 7 + 2] = &perform_HRR_90_72;
      hrrfunc[ANG_HRR_END * 7 + 3] = &perform_HRR_a0_73;
      hrrfunc[ANG_HRR_END * 7 + 4] = &perform_HRR_b0_74;
      hrrfunc[ANG_HRR_END * 7 + 5] = &perform_HRR_c0_75;
      hrrfunc[ANG_HRR_END * 7 + 6] = &perform_HRR_d0_76;
      hrrfunc[ANG_HRR_END * 7 + 7] = &perform_HRR_e0_77;
#endif
}


CHRRList::CHRRList() {
      hrrfunc[ANG_HRR_END * 1 + 1] = &perform_HRR_20_11;
      hrrfunc[ANG_HRR_END * 2 + 1] = &perform_HRR_30_21;
      hrrfunc[ANG_HRR_END * 2 + 2] = &perform_HRR_40_22;
      hrrfunc[ANG_HRR_END * 3 + 1] = &perform_HRR_40_31;
      hrrfunc[ANG_HRR_END * 3 + 2] = &perform_HRR_50_32;
      hrrfunc[ANG_HRR_END * 3 + 3] = &perform_HRR_60_33;
      hrrfunc[ANG_HRR_END * 4 + 1] = &perform_HRR_50_41;
      hrrfunc[ANG_HRR_END * 4 + 2] = &perform_HRR_60_42;
      hrrfunc[ANG_HRR_END * 4 + 3] = &perform_HRR_70_43;
      hrrfunc[ANG_HRR_END * 4 + 4] = &perform_HRR_80_44;
      hrrfunc[ANG_HRR_END * 5 + 1] = &perform_HRR_60_51;
      hrrfunc[ANG_HRR_END * 5 + 2] = &perform_HRR_70_52;
      hrrfunc[ANG_HRR_END * 5 + 3] = &perform_HRR_80_53;
      hrrfunc[ANG_HRR_END * 5 + 4] = &perform_HRR_90_54;
      hrrfunc[ANG_HRR_END * 5 + 5] = &perform_HRR_a0_55;
      hrrfunc[ANG_HRR_END * 6 + 1] = &perform_HRR_70_61;
      hrrfunc[ANG_HRR_END * 6 + 2] = &perform_HRR_80_62;
      hrrfunc[ANG_HRR_END * 6 + 3] = &perform_HRR_90_63;
      hrrfunc[ANG_HRR_END * 6 + 4] = &perform_HRR_a0_64;
      hrrfunc[ANG_HRR_END * 6 + 5] = &perform_HRR_b0_65;
      hrrfunc[ANG_HRR_END * 6 + 6] = &perform_HRR_c0_66;
#ifdef COMPILE_J_ORB
      hrrfunc[ANG_HRR_END * 7 + 1] = &perform_HRR_80_71;
      hrrfunc[ANG_HRR_END * 7 + 2] = &perform_HRR_90_72;
      hrrfunc[ANG_HRR_END * 7 + 3] = &perform_HRR_a0_73;
      hrrfunc[ANG_HRR_END * 7 + 4] = &perform_HRR_b0_74;
      hrrfunc[ANG_HRR_END * 7 + 5] = &perform_HRR_c0_75;
      hrrfunc[ANG_HRR_END * 7 + 6] = &perform_HRR_d0_76;
      hrrfunc[ANG_HRR_END * 7 + 7] = &perform_HRR_e0_77;
#endif
}
