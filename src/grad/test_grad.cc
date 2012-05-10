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
#include <src/grad/goverlapbatch.h>
#include <tuple>
#include <src/grad/gkineticbatch.h>

using namespace std;

void test_grad(shared_ptr<Reference> ref) {

  cout << "  testing grad.." << endl;

  // first construct explicitly grad1e object and store Natom*Nbasis**2 data.
  // This allows us to compare directly the integrals.

  shared_ptr<Geometry> geom = ref->geom();
  const int natom = geom->natom();

  vector<shared_ptr<Matrix1e> > grad1e, grado;
  cout << "  * constructing " << natom*3 << " * " << geom->nbasis() << " * " << geom->nbasis() << endl; 
  for (int i = 0; i != natom*3; ++i) {
    shared_ptr<Matrix1e> a(new Matrix1e(geom));
    grad1e.push_back(a);
    shared_ptr<Matrix1e> b(new Matrix1e(geom));
    grado.push_back(b);
  }

  const vector<shared_ptr<Atom> > atoms = geom->atoms(); 
  const vector<vector<int> > offsets = geom->offsets();
  const int nbasis = geom->nbasis();

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
          const int ang1 = b1->angular_number();
          const int ang0 = b0->angular_number();
          const int s = (ang1+1)*(ang1+2)*(ang0+1)*(ang0+2)/4 * b1->num_primitive()*b0->num_primitive();

          {
            GNAIBatch batch2(input, geom, tie(iatom1, iatom0));
            batch2.compute();
            const double* ndata = batch2.data();
            for (int ia = 0; ia != natom*3; ++ia) {
              for (int i = offset0, cnt = 0; i != dimb0 + offset0; ++i) {
                for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                  grad1e[ia]->data(i*nbasis+j) += ndata[cnt+s*(ia)];
                }
              }
            }
          }
          {
            GKineticBatch batch(input, geom);
            const double* kdata = batch.data();
            batch.compute();
            for (int i = offset0, cnt = 0; i != dimb0 + offset0; ++i) {
              for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                int jatom0 = batch.swap01() ? iatom1 : iatom0;
                int jatom1 = batch.swap01() ? iatom0 : iatom1;
                grad1e[3*jatom1+0]->data(i*nbasis+j) += kdata[cnt];
                grad1e[3*jatom1+1]->data(i*nbasis+j) += kdata[cnt+s];
                grad1e[3*jatom1+2]->data(i*nbasis+j) += kdata[cnt+s*2];
                grad1e[3*jatom0+0]->data(i*nbasis+j) += kdata[cnt+s*3];
                grad1e[3*jatom0+1]->data(i*nbasis+j) += kdata[cnt+s*4];
                grad1e[3*jatom0+2]->data(i*nbasis+j) += kdata[cnt+s*5];
              }
            }
          }
          {
            GOverlapBatch batch(input, geom);
            const double* odata = batch.data();
            batch.compute();
            for (int i = offset0, cnt = 0; i != dimb0 + offset0; ++i) {
              for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                int jatom0 = batch.swap01() ? iatom1 : iatom0;
                int jatom1 = batch.swap01() ? iatom0 : iatom1;
                grado[3*jatom1+0]->data(i*nbasis+j) += odata[cnt];
                grado[3*jatom1+1]->data(i*nbasis+j) += odata[cnt+s];
                grado[3*jatom1+2]->data(i*nbasis+j) += odata[cnt+s*2];
                grado[3*jatom0+0]->data(i*nbasis+j) += odata[cnt+s*3];
                grado[3*jatom0+1]->data(i*nbasis+j) += odata[cnt+s*4];
                grado[3*jatom0+2]->data(i*nbasis+j) += odata[cnt+s*5];
              }
            }
          }

        }
      } 
    }
  }

  shared_ptr<Matrix1e> coeff_occ = ref->coeff()->slice(0,ref->nocc());
  shared_ptr<Matrix1e> rdm1 = ref->rdm1();
  shared_ptr<Matrix1e> erdm1(new Matrix1e(*coeff_occ % *ref->coeff()->form_weighted_density_rhf(ref->nocc(), ref->eig()) * *coeff_occ));

  vector<double> gradv(3*natom);
  for (int i = 0; i != natom*3; ++i) {
    Matrix1e tmp(*coeff_occ % *grad1e[i] * *coeff_occ);
    Matrix1e tmp2(*coeff_occ % *grado[i] * *coeff_occ);
//  gradv[i] = tmp.ddot(rdm1) - tmp2.ddot(erdm1);
    cout << setprecision(10) << tmp.ddot(rdm1) << " " << tmp2.ddot(erdm1) << endl; 
  }
  // the derivative of Vnuc
  auto giter = gradv.begin();
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, giter+=3) {
    const double ax = (*aiter)->position(0);
    const double ay = (*aiter)->position(1);
    const double az = (*aiter)->position(2);
    const double ac = (*aiter)->atom_number();
    for (auto biter = atoms.begin(); biter != atoms.end(); ++biter) {
      if (aiter == biter) continue;
      const double bx = (*biter)->position(0);
      const double by = (*biter)->position(1);
      const double bz = (*biter)->position(2);
      const double c = (*biter)->atom_number() * ac;
      const double dist = sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz));
      *(giter+0) += c*(bx-ax)/(dist*dist*dist);
      *(giter+1) += c*(by-ay)/(dist*dist*dist);
      *(giter+2) += c*(bz-az)/(dist*dist*dist);
    }
  }
  for (int i = 0; i != natom*3; ++i) {
    cout << setprecision(10) << setw(20) << fixed << gradv[i] << endl;
  }


  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref->geom()->df()->compute_half_transform(coeff_occ->data(), ref->nocc());
  // (J^-1)_qp(p|ij)
  shared_ptr<DF_Full> qij =  half->compute_second_transform(coeff_occ->data(), ref->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_closed_2RDM();
  unique_ptr<double[]> qq = qij->form_aux_2index(qijd);

  shared_ptr<DF_AO> qrs = qijd->back_transform(ref->coeff()->data())->back_transform(ref->coeff()->data());

  unique_ptr<double[]> grad2e(new double[3*natom]);
  unique_ptr<double[]> grad2ex(new double[3*natom]);
  fill(grad2e.get(), grad2e.get()+3*natom, 0.0);
  fill(grad2ex.get(), grad2ex.get()+3*natom, 0.0);

  vector<shared_ptr<Atom> > aux_atoms = geom->aux_atoms();
  vector<vector<int> > aux_offsets = geom->aux_offsets();

  // loop over atoms
  for (int iatom0 = 0; iatom0 != geom->natom(); ++iatom0) {
    const shared_ptr<Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();

    for (int iatom1 = 0; iatom1 != geom->natom(); ++iatom1) {
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

          for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
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
                grad2e[3*jatom[i/3]+i%3] += sum;
              }
            }
          }
        }

      }
    }
  }

  for (int iatom0 = 0; iatom0 != geom->natom(); ++iatom0) {
    const shared_ptr<Atom> catom0 = aux_atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = aux_offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();
    for (int iatom1 = 0; iatom1 != geom->natom(); ++iatom1) {
      const shared_ptr<Atom> catom1 = aux_atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = aux_offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();
      // dummy shell
      const shared_ptr<Shell> b3(new Shell(shell1.front()->spherical()));
      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<Shell> b0 = shell0[ibatch0];

        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
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
            grad2ex[3*jatom[i/3]+i%3] -= sum;
          }
        }
      }
    }
  }

  cout << endl;
  for (int i = 0; i != natom*3; ++i) {
    cout << setw(3) << i << setprecision(10) << setw(20) << fixed << grad2e[i] << setw(20) << grad2ex[i] << endl;
  }

}
