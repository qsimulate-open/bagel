//
// BAGEL - Parallel electron correlation program.
// Filename: pfock.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/pscf/pfock.h>
#include <src/util/f77.h>
#include <cstring>

using namespace std;
using namespace bagel;

typedef shared_ptr<PFock> RefPFock;
typedef shared_ptr<PHcore> RefPHcore;
typedef shared_ptr<PGeometry> RefPGeometry;
typedef shared_ptr<PMatrix1e> RefPMatrix1e;

PFock::PFock(const RefPGeometry g, const RefPFock prev, const RefPMatrix1e den, const vector<double>& shw, const int th, const bool dir)
 : PMatrix1e(g), previous_(prev), density_(den), schwarz_(shw), S2_(th), direct_(dir) {

  const int unit = 1;
  const complex<double> one(1.0, 0.0);

  // First construct the 2-e Fock matrices
  pfock_two_electron_part();

  PMatrix1e tmp = this->ft();
  zcopy_(&totalsize_, tmp.data()->front(), &unit, data_->front(), &unit);

  // Finally, add the old Fock matrices
  zaxpy_(&totalsize_, &one, prev->data()->front(), &unit, data_->front(), &unit);


}

PFock::PFock(const RefPGeometry g, const RefPFock prev, const RefPMatrix1e den, const vector<double>& shw, const int th, const bool dir, shared_ptr<PCompFile<ERIBatch>> fl)
 : PMatrix1e(g), previous_(prev), density_(den), schwarz_(shw), S2_(th), direct_(dir) {
  assert(!direct());

  file_ = fl;

  const int unit = 1;
  const complex<double> one(1.0, 0.0);
  // First construct the 2-e Fock matrices
  pfock_two_electron_part();
  PMatrix1e tmp = this->ft();
  zcopy_(&totalsize_, tmp.data()->front(), &unit, data_->front(), &unit);
  // Finally, add the old Fock matrices
  zaxpy_(&totalsize_, &one, prev->data()->front(), &unit, data_->front(), &unit);
}

PFock::PFock(const RefPGeometry g, const RefPHcore hcore) : PMatrix1e(g) {

  PMatrix1e hcore_m = hcore->ft();
  const int unit = 1;
  zcopy_(&totalsize_, hcore_m.data()->front(), &unit, data_->front(), &unit);

}


PFock::~PFock() {

}


