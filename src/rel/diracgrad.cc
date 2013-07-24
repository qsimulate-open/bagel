//
// BAGEL - Parallel electron correlation program.
// Filename: diracgrad.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/rel/alpha.h>
#include <src/rel/dfock.h>
#include <src/rel/reldffull.h>
#include <src/rel/cdmatrix.h>
#include <src/integral/rys/gsmallnaibatch.h>

using namespace std;
using namespace bagel;


template<>
shared_ptr<GradFile> GradEval<Dirac>::compute() {
  geom_ = task_->geom();

  Timer timer;
  // density matrix
  shared_ptr<const RelReference> ref = dynamic_pointer_cast<const RelReference>(ref_);
  shared_ptr<const ZMatrix> coeff = ref->relcoeff()->slice(0, ref->nocc());
  auto den = make_shared<const ZMatrix>(*coeff ^ *coeff);

  // energy-weighted density matrix
  shared_ptr<ZMatrix> ecoeff = coeff->copy();
  const vector<double>& eig = ref->eig();
  for (int i = 0; i != ref->nocc(); ++i)
    zscal_(ecoeff->ndim(), eig[i], ecoeff->element_ptr(0, i), 1);
  auto eden = make_shared<const ZMatrix>(*coeff ^ *ecoeff);

  const int nbasis = geom_->nbasis();

  // NAI density (L+L+, L-L-)
  shared_ptr<ZMatrix> nden  =  den->get_submatrix(0, 0, nbasis, nbasis);
                     *nden += *den->get_submatrix(nbasis, nbasis, nbasis, nbasis);
  // kinetic density [den] 2*(S+L+, S-L-) - (S+S+, S-S-); [eden] -(S+S+, S-S-)/2c^2
  shared_ptr<ZMatrix> kden  =  den->get_submatrix(0, 2*nbasis, nbasis, nbasis);
                     *kden += *den->get_submatrix(nbasis, 3*nbasis, nbasis, nbasis);
                     *kden *= complex<double>(2.0);
                     *kden -= *den->get_submatrix(2*nbasis, 2*nbasis, nbasis, nbasis);
                     *kden -= *den->get_submatrix(3*nbasis, 3*nbasis, nbasis, nbasis);
  shared_ptr<ZMatrix> lden  = eden->get_submatrix(2*nbasis, 2*nbasis, nbasis, nbasis);
                     *lden +=*eden->get_submatrix(3*nbasis, 3*nbasis, nbasis, nbasis);
                     *lden /= complex<double>(2.0*pow(c__,2));
                     *kden -= *lden;
  // overlap density
  shared_ptr<ZMatrix> sden  = eden->get_submatrix(0, 0, nbasis, nbasis);
                     *sden +=*eden->get_submatrix(nbasis, nbasis, nbasis, nbasis);

  // nden, kden, sden (minus sign is taken care of inside)
  vector<GradTask> task = contract_grad1e(nden->get_real_part(), kden->get_real_part(), sden->get_real_part());

  // small NAI part..
  map<int, shared_ptr<Sigma>> sigma;
  sigma.insert(make_pair(0, make_shared<Sigma>(Comp::X)));
  sigma.insert(make_pair(1, make_shared<Sigma>(Comp::Y)));
  sigma.insert(make_pair(2, make_shared<Sigma>(Comp::Z)));
  auto sp = make_shared<ZMatrix>(4,1,true); sp->element(2,0) = 1;
  auto sm = make_shared<ZMatrix>(4,1,true); sm->element(3,0) = 1;
  map<int, shared_ptr<ZMatrix>> XY{ make_pair(2, sp), make_pair(3, sm) };

  // target data area
  vector<int> xyz{Comp::X, Comp::Y, Comp::Z};
  map<pair<int,int>, shared_ptr<ZMatrix>> mat;
  for (auto& i : xyz)
    for (auto& j : xyz)
      if (i <= j)
        mat.insert(make_pair(make_pair(i,j), make_shared<ZMatrix>(nbasis, nbasis)));

  for (auto& s0 : XY) { // bra
    for (auto& s1 : XY) { // ket
      shared_ptr<ZMatrix> data = den->get_submatrix(s0.first*nbasis, s1.first*nbasis, nbasis, nbasis);
      for (auto& w0 : sigma) {
        for (auto& w1 : sigma) {
          auto tmp = make_shared<ZMatrix>((*w0.second * *s0.second) % (*w1.second * *s1.second));
          const complex<double> c = tmp->element(0,0);
          const int small = min(w0.first, w1.first);
          const int large = max(w0.first, w1.first);
          mat[make_pair(small, large)]->zaxpy(c, data);
        }
      }
    }
  }
  assert(mat.size() == 6);
  array<shared_ptr<const Matrix>,6> rmat;
  auto riter = rmat.begin();
  for (auto& i : mat) *riter++ = i.second->get_real_part();
  vector<GradTask> tmp = contract_gradsmall1e(rmat);
  task.insert(task.end(), tmp.begin(), tmp.end());

  // two-electron contributions.
  {
    // make blocks of coefficients
    array<shared_ptr<const Matrix>, 4> rocoeff;
    array<shared_ptr<const Matrix>, 4> iocoeff;
    array<shared_ptr<const Matrix>, 4> trocoeff;
    array<shared_ptr<const Matrix>, 4> tiocoeff;
    for (int i = 0; i != 4; ++i) {
      shared_ptr<const ZMatrix> ocoeff = coeff->get_submatrix(i*geom_->nbasis(), 0, geom_->nbasis(), ref->nocc());
      rocoeff[i] = ocoeff->get_real_part();
      iocoeff[i] = ocoeff->get_imag_part();
      trocoeff[i] = rocoeff[i]->transpose();
      tiocoeff[i] = iocoeff[i]->transpose();
    }
    // (0) get AO integras
    // get individual df dist objects for each block and add df to dfs
    vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());

    // (1) make RelDF objects from AO integrals
    list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

    // (2) first-transform
    list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, rocoeff, iocoeff);

    // (3) split and factorize
    list<shared_ptr<RelDFHalf>> half_complex_exch;
    for (auto& i : half_complex) {
      list<shared_ptr<RelDFHalf>> tmp = i->split(false);
      half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
    }
    half_complex.clear();
    DFock::factorize(half_complex_exch);

    // (4) compute C matrix
    shared_ptr<CDMatrix> cd;
    for (auto& j : half_complex_exch) {
      for (auto& i : j->basis()) {
        if (cd) {
          *cd += CDMatrix(j, i, trocoeff, tiocoeff, geom_->df()->data2(), false /* J^-1 multiplied */);
        } else {
          cd = make_shared<CDMatrix>(j, i, trocoeff, tiocoeff, geom_->df()->data2(), false /* J^-1 multiplied */);
        }
      }
    }
    cd->localize();

    // (5) compute (gamma|ij)
    list<shared_ptr<RelDFFull>> dffull;
    for (auto& i : half_complex_exch)
      dffull.push_back(make_shared<RelDFFull>(i, rocoeff, iocoeff));
    DFock::factorize(dffull);
    assert(dffull.size() == 1);

    // (6) two-index gamma
    shared_ptr<Matrix> cdr = cd->get_real_part(); 
    shared_ptr<Matrix> cdi = cd->get_imag_part(); 
    auto gamma2 = make_shared<Matrix>((*cdr ^ *cdr) - (*cdi ^ *cdi) - *dffull.front()->form_aux_2index_real());
    gamma2->print();
  }


  // compute
  TaskQueue<GradTask> tq(task);
  tq.compute(resources__->max_num_threads());

  // allreduce
  mpi__->allreduce(grad_->data()->data(), grad_->size());

  // adds nuclear contributions
  *grad_->data() += *geom_->compute_grad_vnuc();

  grad_->print();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad_;
}

