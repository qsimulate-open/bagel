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

#include <src/grad/gradeval_base.h>

using namespace std;

void test_grad(shared_ptr<Reference> ref) {

  cout << "  testing grad.." << endl;

  const shared_ptr<const Geometry> geom = ref->geom();
  const int natom = geom->natom();
  // target quantity here ... ==========
  size_t start = ::clock(); 
  vector<double> grad(3*natom);

  shared_ptr<Grad1eFile> grad1e, grado;
  GradEval_base g1(geom);
  tie(grad1e, grado) = g1.data();

  shared_ptr<Matrix1e> coeff_occ = ref->coeff()->slice(0,ref->nocc());
  shared_ptr<Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref->rdm1() ^ *coeff_occ));
  shared_ptr<Matrix1e> erdm1 = ref->coeff()->form_weighted_density_rhf(ref->nocc(), ref->eig());

  vector<double> gvnuc = geom->compute_grad_vnuc();
  auto giter = gvnuc.begin();
  for (int i = 0; i != natom*3; ++i, ++giter)
    grad[i] = (*grad1e)[i]->ddot(rdm1) - (*grado)[i]->ddot(erdm1) + *giter;

  cout << left << setw(50) << "  * One body part computed with " << setprecision(3) << right << setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  start = ::clock();


  const vector<shared_ptr<Atom> > atoms = geom->atoms(); 
  const vector<vector<int> > offsets = geom->offsets();
  const int nbasis = geom->nbasis();


  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref->geom()->df()->compute_half_transform(coeff_occ->data(), ref->nocc());
  // (J^-1)_qp(p|ij)
  shared_ptr<DF_Full> qij =  half->compute_second_transform(coeff_occ->data(), ref->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_closed_2RDM();
  unique_ptr<double[]> qq = qij->form_aux_2index(qijd);

  shared_ptr<DF_AO> qrs = qijd->back_transform(ref->coeff()->data())->back_transform(ref->coeff()->data());

  cout << left << setw(50) << "  * Two body part (trans) computed with " << setprecision(3) << right << setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  start = ::clock();


  vector<shared_ptr<Atom> > aux_atoms = geom->aux_atoms();
  vector<vector<int> > aux_offsets = geom->aux_offsets();

  // loop over atoms (using symmetry b0 <-> b1)
  for (int iatom0 = 0; iatom0 != geom->natom(); ++iatom0) {
    const shared_ptr<Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();

    for (int iatom1 = iatom0; iatom1 != geom->natom(); ++iatom1) {
      const shared_ptr<Atom> catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();

      for (int iatom2 = 0; iatom2 != geom->natom(); ++iatom2) {
        const shared_ptr<Atom> catom2 = aux_atoms[iatom2];
        const int numshell2 = catom2->shells().size();
        const vector<int> coffset2 = aux_offsets[iatom2];
        const vector<shared_ptr<Shell> > shell2 = catom2->shells();


        // dummy shell
        const shared_ptr<Shell> b3(new Shell(shell2.front()->spherical()));

        for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
          const int offset0 = coffset0[ibatch0]; 
          shared_ptr<Shell> b0 = shell0[ibatch0];

          for (int ibatch1 = (iatom0 == iatom1 ? ibatch0 : 0); ibatch1 != numshell1; ++ibatch1) {
            const int offset1 = coffset1[ibatch1]; 
            shared_ptr<Shell> b1 = shell1[ibatch1];

            for (int ibatch2 = 0; ibatch2 != numshell2; ++ibatch2) {
              const int offset2 = coffset2[ibatch2]; 
              shared_ptr<Shell> b2 = shell2[ibatch2];

              vector<shared_ptr<Shell> > input;
              input.push_back(b3);
              input.push_back(b2);
              input.push_back(b1);
              input.push_back(b0);

              // pointer to stack
              GradBatch gradbatch(input, 0.0);
              gradbatch.compute();
              const size_t block = gradbatch.size_block();

              // unfortunately the convention is different...
              int jatom[4] = {-1, iatom2, iatom1, iatom0};
              if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
              if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
              if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

              for (int i = 0; i != 12; ++i) {
                // if this is a dummy atom
                if (jatom[i/3] < 0) continue;

                const double* ppt = gradbatch.data() + i*block;
                double sum = 0.0;
                for (int j0 = offset0; j0 != offset0 + b0->nbasis(); ++j0) {  
                  for (int j1 = offset1; j1 != offset1 + b1->nbasis(); ++j1) {  
                    for (int j2 = offset2; j2 != offset2 + b2->nbasis(); ++j2, ++ppt) {  
                      sum += *ppt * *qrs->ptr(j2, j1, j0);
                    }
                  }
                }
                grad[3*jatom[i/3]+i%3] += sum * (iatom0 == iatom1 && ibatch0 == ibatch1 ? 1.0 : 2.0);
              }
            }
          }
        }

      }
    }
  }

  cout << left << setw(50) << "  * Two body part (3-center) computed with " << setprecision(3) << right << setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  start = ::clock();

  // using symmetry (b0 <-> b1)
  for (int iatom0 = 0; iatom0 != geom->natom(); ++iatom0) {
    const shared_ptr<Atom> catom0 = aux_atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = aux_offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();
    for (int iatom1 = iatom0; iatom1 != geom->natom(); ++iatom1) {
      const shared_ptr<Atom> catom1 = aux_atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = aux_offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();
      // dummy shell
      const shared_ptr<Shell> b3(new Shell(shell1.front()->spherical()));
      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<Shell> b0 = shell0[ibatch0];

        for (int ibatch1 = (iatom0 == iatom1 ? ibatch0 : 0); ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          shared_ptr<Shell> b1 = shell1[ibatch1];

          vector<shared_ptr<Shell> > input;
          input.push_back(b1);
          input.push_back(b3);
          input.push_back(b0);
          input.push_back(b3);

          // pointer to stack
          GradBatch gradbatch(input, 0.0);
          gradbatch.compute();
          const size_t block = gradbatch.size_block();

          // unfortunately the convention is different...
          int jatom[4] = {iatom1, -1, iatom0, -1};
          if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
          if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
          if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

          for (int i = 0; i != 12; ++i) {
            // if this is a dummy atom
            if (jatom[i/3] < 0) continue;

            const double* ppt = gradbatch.data() + i*block;
            double sum = 0.0;
            for (int j0 = offset0; j0 != offset0 + b0->nbasis(); ++j0) {  
              for (int j1 = offset1; j1 != offset1 + b1->nbasis(); ++j1, ++ppt) {  
                sum += *ppt * qq[j1+ref->geom()->naux()*j0];
              }
            }
            grad[3*jatom[i/3]+i%3] -= sum * 0.5 * (iatom0 == iatom1 && ibatch0 == ibatch1 ? 1.0 : 2.0);
          }
        }
      }
    }
  }

  cout << endl << "  * Nuclear energy gradient" << endl << endl;
  for (int i = 0; i != natom; ++i) {
    cout << "    o Atom " << setw(3) << i << endl;
    cout << "        x  " << setprecision(10) << setw(20) << fixed << grad[3*i+0] << endl;
    cout << "        y  " << setprecision(10) << setw(20) << fixed << grad[3*i+1] << endl;
    cout << "        z  " << setprecision(10) << setw(20) << fixed << grad[3*i+2] << endl;
  }

  cout << setw(50) << left << "  * Two body part (2-center) computed with " << setprecision(3) << right << setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  cout << endl;

}


