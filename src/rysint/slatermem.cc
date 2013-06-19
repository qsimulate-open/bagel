//
// BAGEL - Parallel electron correlation program.
// Filename: slatermem.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include <src/rysint/slatermem.h>

using namespace std;
using namespace bagel;

// will be used by root*.cc
namespace bagel {
  const SlaterMem slatermem__;
}

SlaterMem::SlaterMem() {
  for (int i = 1; i != 14; ++i) {
    datax_.push_back(std::unique_ptr<double[]>(new double[19600*i]));
    dataw_.push_back(std::unique_ptr<double[]>(new double[19600*i]));
  }
  fill1();
  fill2();
  fill3();
  fill4();
  fill5();
  fill6();
  fill7();
  fill8();
  fill9();
  fill10();
  fill11();
  fill12();
  fill13();
}
