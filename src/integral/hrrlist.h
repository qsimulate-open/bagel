//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hrrlist.h
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

//
// stand alone class for HRR functions.
// can be called from outside through hrrfunc_call.
//

#ifndef __SRC_RYSINT_HRRLIST_H
#define __SRC_RYSINT_HRRLIST_H

#include <complex>
#include <functional>
#include <array>
#include <src/util/constants.h>

namespace bagel {

struct HRRList {
  HRRList();

  static void perform_HRR_20_11(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_30_21(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_40_22(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_40_31(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_50_32(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_50_41(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_60_33(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_60_42(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_60_51(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_70_43(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_70_52(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_70_61(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_80_44(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_80_53(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_80_62(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_90_54(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_90_63(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_a0_55(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_a0_64(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_b0_65(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_c0_66(const int, const double*, const std::array<double,3>&, double*);
#ifdef COMPILE_J_ORB
  static void perform_HRR_80_71(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_90_72(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_a0_73(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_b0_74(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_c0_75(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_d0_76(const int, const double*, const std::array<double,3>&, double*);
  static void perform_HRR_e0_77(const int, const double*, const std::array<double,3>&, double*);
#endif

  void hrrfunc_call(const unsigned int i, const int a0, const double* a1, const std::array<double,3>& a2, double* a3) const {
    return hrrfunc[i](a0, a1, a2, a3);
  }

  std::function<void (const int, const double*, const std::array<double,3>&, double*)> hrrfunc[ANG_HRR_END * ANG_HRR_END];
};


struct CHRRList {
  CHRRList();

  static void perform_HRR_20_11(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_30_21(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_40_22(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_40_31(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_50_32(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_50_41(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_60_33(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_60_42(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_60_51(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_70_43(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_70_52(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_70_61(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_80_44(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_80_53(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_80_62(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_90_54(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_90_63(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_a0_55(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_a0_64(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_b0_65(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_c0_66(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
#ifdef COMPILE_J_ORB
  static void perform_HRR_80_71(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_90_72(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_a0_73(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_b0_74(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_c0_75(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_d0_76(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
  static void perform_HRR_e0_77(const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*);
#endif

  void hrrfunc_call(const unsigned int i, const int a0, const std::complex<double>* a1, const std::array<double,3>& a2, std::complex<double>* a3) const {
    return hrrfunc[i](a0, a1, a2, a3);
  }

  std::function<void (const int, const std::complex<double>*, const std::array<double,3>&, std::complex<double>*)> hrrfunc[ANG_HRR_END * ANG_HRR_END];
};

}

#endif
