//
// Newint - Parallel electron correlation program.
// Filename: momentum.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/prop/momentum.h>
#include <src/osint/momentbatch.h>
#include <iomanip>

using namespace std;

Momentum::Momentum(shared_ptr<const Geometry> g) : geom_(g) {

}


Momentum::~Momentum() {

}


array<shared_ptr<Matrix1e>, 3> Momentum::compute() const {
  const int natom = geom_->natom();

  const vector<shared_ptr<Atom> > atoms = geom_->atoms(); 
  const vector<vector<int> > offsets = geom_->offsets();
  const int nbasis = geom_->nbasis();

  const shared_ptr<Matrix1e> mat0(new Matrix1e(geom_));
  const shared_ptr<Matrix1e> mat1(new Matrix1e(geom_));
  const shared_ptr<Matrix1e> mat2(new Matrix1e(geom_));

  // TODO perhaps we could reduce operation by a factor of 2
  for (int iatom0 = 0; iatom0 != natom; ++iatom0) {
    const shared_ptr<Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<const Shell> > shell0 = catom0->shells();

    for (int iatom1 = 0; iatom1 != natom; ++iatom1) {
      const shared_ptr<Atom> catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<shared_ptr<const Shell> > shell1 = catom1->shells();

      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<const Shell> b0 = shell0[ibatch0];
        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          shared_ptr<const Shell> b1 = shell1[ibatch1];

          vector<shared_ptr<const Shell> > input = {{b1, b0}};
          MomentBatch mom(input);
          mom.compute();

          const int dimb1 = input[0]->nbasis(); 
          const int dimb0 = input[1]->nbasis(); 
          const double* dat0 = mom.data();
          const double* dat1 = mom.data() + mom.size_block();
          const double* dat2 = mom.data() + mom.size_block()*2;
          for (int i = offset0; i != dimb0 + offset0; ++i) {
            for (int j = offset1; j != dimb1 + offset1; ++j, ++dat0, ++dat1, ++dat2) {
              mat0->element(j,i) = *dat0;
              mat1->element(j,i) = *dat1;
              mat2->element(j,i) = *dat2;
            }
          }

        }
      }
    }
  }

  return array<shared_ptr<Matrix1e>,3>{{mat0, mat1, mat2}};
}
