//
// Newint - Parallel electron correlation program.
// Filename: carsphlist.h
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

  
#ifndef __rysint_carsphlist_h
#define __rysint_carsphlist_h

#include <src/rysint/macros.h>
 
struct CarSphList {
  CarSphList();
  ~CarSphList();
     
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

  void carsphfunc_call(const unsigned int i, const int a0, const double* a1, double* a2) {
    return (carsphfunc[i])(a0, a1, a2);
  };

  void (*carsphfunc[ANG_HRR_END * ANG_HRR_END])(const int, const double*, double*);
};

#endif
