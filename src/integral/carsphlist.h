//
// BAGEL - Parallel electron correlation program.
// Filename: carsphlist.h
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


#ifndef __SRC_RYSINT_CARSPHLIST_H
#define __SRC_RYSINT_CARSPHLIST_H

#include <complex>
#include <functional>
#include <src/util/constants.h>

namespace bagel {

struct CarSphList {
  CarSphList();

  static void carsph_00(const int, const double*, double*);
  static void carsph_10(const int, const double*, double*);
  static void carsph_20(const int, const double*, double*);
  static void carsph_30(const int, const double*, double*);
  static void carsph_40(const int, const double*, double*);
  static void carsph_50(const int, const double*, double*);
  static void carsph_60(const int, const double*, double*);
  static void carsph_11(const int, const double*, double*);
  static void carsph_21(const int, const double*, double*);
  static void carsph_31(const int, const double*, double*);
  static void carsph_41(const int, const double*, double*);
  static void carsph_51(const int, const double*, double*);
  static void carsph_61(const int, const double*, double*);
  static void carsph_22(const int, const double*, double*);
  static void carsph_32(const int, const double*, double*);
  static void carsph_42(const int, const double*, double*);
  static void carsph_52(const int, const double*, double*);
  static void carsph_62(const int, const double*, double*);
  static void carsph_33(const int, const double*, double*);
  static void carsph_43(const int, const double*, double*);
  static void carsph_53(const int, const double*, double*);
  static void carsph_63(const int, const double*, double*);
  static void carsph_44(const int, const double*, double*);
  static void carsph_54(const int, const double*, double*);
  static void carsph_64(const int, const double*, double*);
  static void carsph_55(const int, const double*, double*);
  static void carsph_65(const int, const double*, double*);
  static void carsph_66(const int, const double*, double*);

  void carsphfunc_call(const unsigned int i, const int a0, const double* a1, double* a2) const {
    return (carsphfunc[i])(a0, a1, a2);
  }

  std::function<void (const int, const double*, double*)> carsphfunc[ANG_HRR_END * ANG_HRR_END];
};


struct CCarSphList {
  CCarSphList();

  static void carsph_00(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_10(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_20(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_30(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_40(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_50(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_60(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_11(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_21(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_31(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_41(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_51(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_61(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_22(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_32(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_42(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_52(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_62(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_33(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_43(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_53(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_63(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_44(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_54(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_64(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_55(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_65(const int, const std::complex<double>*, std::complex<double>*);
  static void carsph_66(const int, const std::complex<double>*, std::complex<double>*);

  void carsphfunc_call(const unsigned int i, const int a0, const std::complex<double>* a1, std::complex<double>* a2) const {
    return (carsphfunc[i])(a0, a1, a2);
  }

  std::function<void (const int, const std::complex<double>*, std::complex<double>*)> carsphfunc[ANG_HRR_END * ANG_HRR_END];
};

}

#endif
