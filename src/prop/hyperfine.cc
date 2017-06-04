//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hyperfine.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/prop/hyperfine.h>
#include <src/mat1e/fermicontact.h>
#include <src/mat1e/spindipole.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap;

HyperFine::HyperFine(shared_ptr<const Geometry> geom, shared_ptr<const Matrix> den, const int s, const string jobname, vector<int> select)
 : geom_(geom), den_(den), jobname_(jobname), select_(select), s_(s) {

}


void HyperFine::compute() const {
  const string indent = "      ";
  cout << "    * Hyperfine coupling constants (" << jobname_ << ")" << endl;

  vector<shared_ptr<const Atom>> atoms = geom_->atoms();
  int cnt = 0;
  for (auto& i : atoms) {
    if (!select_.empty() && find(select_.begin(), select_.end(), cnt) == select_.end()) continue;

    if (atommap.hfcc_exists(i->name())) {
      cout << indent << "Atom: " << setw(4) << cnt << endl;
      cout << indent << "  Spin dipole" << setprecision(10) << endl;
      SpinDipole mat(geom_, i, s_);
      for (int j = 0; j != 6; ++j)
        cout << setw(22) << mat.data(j)->dot_product(den_) << endl;

      cout << indent << "  Fermi contact" << endl;
      FermiContact mat2(geom_, i, s_);
      cout << setw(22) << mat2.dot_product(den_) << endl << endl;
    }

    ++cnt;
  }
}
