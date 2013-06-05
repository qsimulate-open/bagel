//
// BAGEL - Parallel electron correlation program.
// Filename: process.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <fstream>
#include <sstream>
#include <src/parallel/process.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

Process::Process() : print_level_(3) {
  if (mpi__->rank() != 0) {
    cout_orig = cout.rdbuf();
    cout.rdbuf(ss_.rdbuf());
  }
}

Process::~Process() {
  if (mpi__->rank() != 0)
    cout.rdbuf(cout_orig);
}


void Process::cout_on()  const { if (mpi__->rank() != 0) cout.rdbuf(cout_orig); }
void Process::cout_off() const { if (mpi__->rank() != 0) cout.rdbuf(ss_.rdbuf()); }
