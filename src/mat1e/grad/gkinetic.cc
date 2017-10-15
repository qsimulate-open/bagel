//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gkinetic.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#include <src/mat1e/grad/gkinetic.h>
#include <src/integral/os/gkineticbatch.h>

using namespace std;
using namespace bagel;

GKinetic::GKinetic(shared_ptr<const Molecule> mol) : GMatrix1e(mol) {

	init(mol);
	// fill_upper();

}


void GKinetic::computebatch(const array<shared_ptr<const Shell>,2>& s, const vector<int>& a, const vector<int>& o, shared_ptr<const Molecule> m) {

	// s = [b1, b0]
	assert(s.size() == 2);

	GKineticBatch batch(s);
	batch.compute();

	const int dimb1 = s[0]->nbasis(), dimb0 = s[1]->nbasis();
	const int iatom0 = a[0], iatom1 = a[1], natom = m->natom();
	const int offsetb0 = o[0], offsetb1 = o[1];
	const int jatom0 = batch.swap01() ? iatom1 : iatom0, jatom1 = batch.swap01() ? iatom0 : iatom1;

	for (int k = 0; k != 3; ++k) {
		for (int i = offsetb0, j = 0; i != offsetb0 + dimb0 ; ++i, ++j) {
			blas::ax_plus_y_n(1, batch.data() + k * batch.size_block() + j * dimb1, dimb1, matrices_[jatom1]->element_ptr(offsetb1, i));
			blas::ax_plus_y_n(1, batch.data() + (k + 3) * batch.size_block() + j * dimb1, dimb1, matrices_[jatom0]->element_ptr(offsetb1, i));
		}
	}

}


