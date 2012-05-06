//
// Newint - Parallel electron correlation program.
// Filename: test_grad.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/df/fit.h>
#include <src/wfn/reference.h>
#include <iostream>
#include <src/osint/kineticbatch.h>
#include <src/grad/gnaibatch.h>
#include <tuple>

using namespace std;

void test_grad(shared_ptr<Reference> ref) {

  cout << "  testing grad.." << endl;
#if 0
  shared_ptr<Geometry> g = ref->geom();
  shared_ptr<DensityFit> grad(dynamic_cast<DensityFit*>(new GradFit(g->nbasis(), g->naux(), g->atoms(), g->offsets(),
                                                                    g->aux_atoms(), g->aux_offsets(), 0.0, false)));
#endif

  // first construct explicitly grad1e object and store Natom*Nbasis**2 data.
  // This allows us to compare directly the integrals.

  shared_ptr<Geometry> geom_ = ref->geom();
  const int natom = geom_->natom();

  vector<shared_ptr<Matrix1e> > grad1e;
  cout << "  * constructing " << natom*3 << " * " << geom_->nbasis() << " * " << geom_->nbasis() << endl; 
  for (int i = 0; i != natom*3; ++i) {
    shared_ptr<Matrix1e> a(new Matrix1e(geom_));
    grad1e.push_back(a);
  }

  const vector<shared_ptr<Atom> > atoms = geom_->atoms(); 
  const vector<vector<int> > offsets = geom_->offsets();
  const int nbasis = geom_->nbasis();

  // only lower half will be stored
  for (int iatom0 = 0; iatom0 != natom; ++iatom0) {
    // iatom1 = iatom1;
    const shared_ptr<Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();

    for (int iatom1 = 0; iatom1 != natom; ++iatom1) {
      const shared_ptr<Atom> catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();

      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<Shell> b0 = shell0[ibatch0];
        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          shared_ptr<Shell> b1 = shell1[ibatch1];
          vector<shared_ptr<Shell> > input;
          input.push_back(b1);
          input.push_back(b0);

          const int dimb1 = input[0]->nbasis(); 
          const int dimb0 = input[1]->nbasis(); 
          KineticBatch kinetic(input);
          kinetic.compute();
          const double* kdata = kinetic.data();
          GNAIBatch nai(input, geom_, tie(iatom0, iatom1));
          nai.compute();
          const double* ndata = nai.data();

          const int ang1 = b1->angular_number();
          const int ang0 = b0->angular_number();
          const int s = (ang1+1)*(ang1+2)*(ang0+1)*(ang0+2)/4 * b1->num_primitive()*b0->num_primitive();

          for (int ia = 0; ia != natom*3; ++ia) {
            int cnt = 0;
            for (int i = offset0; i != dimb0 + offset0; ++i) {
              for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                grad1e[ia]->data(i*nbasis+j) += ndata[cnt+s*(ia)];
              }
            }
          }

        }
      } 
    }
  } 
  for (int i = 0; i != natom*3; ++i) {
    grad1e[i]->print("", 12);
  }

}
