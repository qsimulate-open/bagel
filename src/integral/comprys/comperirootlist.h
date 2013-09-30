//
// BAGEL - Parallel electron correlation program.
// Filename: comperirootlist.h
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

#ifndef __SRC_COMPRYSINT_ERIROOTLIST_H
#define __SRC_COMPRYSINT_ERIROOTLIST_H

#include <functional>
#include <complex>
#include <src/util/constants.h>

namespace bagel {

  struct ComplexERIRootList  {
    private:
      std::function<void (const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int)> rfunc[RYS_MAX + 1];

      static void complex_eriroot1(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot2(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot3(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot4(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot5(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot6(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot7(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot8(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot9(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot10(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot11(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot12(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);
      static void complex_eriroot13(const std::complex<double>*, std::complex<double>*, std::complex<double>*, const int);


    public:
      ComplexERIRootList() {
        rfunc[1] = &complex_eriroot1;
        rfunc[2] = &complex_eriroot2;
        rfunc[3] = &complex_eriroot3;
        rfunc[4] = &complex_eriroot4;
        rfunc[5] = &complex_eriroot5;
        rfunc[6] = &complex_eriroot6;
        rfunc[7] = &complex_eriroot7;
        rfunc[8] = &complex_eriroot8;
        rfunc[9] = &complex_eriroot9;
        rfunc[10] = &complex_eriroot10;
        rfunc[11] = &complex_eriroot11;
        rfunc[12] = &complex_eriroot12;
        rfunc[13] = &complex_eriroot13;
      }

      void root(const int i, const std::complex<double>* a1, std::complex<double>* a2, std::complex<double>* a3, const int a4) const { rfunc[i](a1, a2, a3, a4); }

  };
 
}

#endif
