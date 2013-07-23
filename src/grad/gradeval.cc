//
// BAGEL - Parallel electron correlation program.
// Filename: gradeval.cc
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


#include <src/integral/rys/gsmallnaibatch.h>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/rel/alpha.h>

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
          auto tmp = make_shared<ZMatrix>((*w0.second->data() * *s0.second) % (*w1.second->data() * *s1.second));
          const complex<double> c = tmp->element(0,0);
          const int small = min(w0.first, w1.first);
          const int large = max(w0.first, w1.first);
          mat[make_pair(small, large)]->zaxpy(c, data);
        }
      }
    }
  }
  array<shared_ptr<const Matrix>,6> dmat;
  auto iter = mat.begin();
  for (auto& i : dmat) {
    i = iter->second->get_real_part();
    ++iter;
  }

  {
    int mpicnt = 0;
    int iatom0 = 0;
    auto oa0 = geom_->offsets().begin();
    for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++oa0, ++iatom0) {
      int iatom1 = 0;
      auto oa1 = geom_->offsets().begin();
      for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++oa1, ++iatom1) {

        auto o0 = oa0->begin();
        for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
          auto o1 = oa1->begin();
          for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

            // static distribution since this is cheap
            if (mpicnt++ % mpi__->size() != mpi__->rank()) continue;

            // TODO change the following to tasks
            array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
            vector<int> atom = {iatom0, iatom1};
            vector<int> offset_ = {*o0, *o1};
            // make six Cartesian blocks
            GSmallNAIBatch batch(input, geom_, tie(iatom1, iatom0));
            batch.compute();
            *grad_ += *batch.compute_gradient(dmat);
#if 0
            const double* ndata = batch.data();
            const int dimb1 = input[0]->nbasis();
            const int dimb0 = input[1]->nbasis();
            const size_t s = batch.size_block();
            const size_t cartblock = s*geom_->natom()*3;
            for (int ia = 0; ia != geom_->natom()*3; ++ia) {
              for (int i = offset_[0], cnt = 0; i != dimb0 + offset_[0]; ++i) {
                for (int j = offset_[1]; j != dimb1 + offset_[1]; ++j, ++cnt) {
                  grad_->data(ia) += ndata[cnt+s*ia            ] * dmat[0]->element(j,i);
                  grad_->data(ia) += ndata[cnt+s*ia+cartblock  ] * dmat[1]->element(j,i);
                  grad_->data(ia) += ndata[cnt+s*ia+cartblock*2] * dmat[2]->element(j,i);
                  grad_->data(ia) += ndata[cnt+s*ia+cartblock*3] * dmat[3]->element(j,i);
                  grad_->data(ia) += ndata[cnt+s*ia+cartblock*4] * dmat[4]->element(j,i);
                  grad_->data(ia) += ndata[cnt+s*ia+cartblock*5] * dmat[5]->element(j,i);
                }
              }
            }
#endif
          }
        }
      }
    }
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

template<>
shared_ptr<GradFile> GradEval<SCF>::compute() {
  assert(task_->dodf());
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_closed_2RDM();
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<UHF>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();
  assert(erdm1 != nullptr);

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<ROHF>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();
  assert(erdm1 != nullptr);

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<KS>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  // ... exchange needs to be scaled.
  shared_ptr<const DFFullDist> qijd = qij->apply_closed_2RDM(task_->func()->scale_ex());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  //- Exchange-correlation part -//
  shared_ptr<const GradFile> ggrad = task_->grid()->compute_xcgrad(task_->func(), coeff_occ);
  *grad += *ggrad;

  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<WernerKnowles>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}

template<>
shared_ptr<GradFile> GradEval<SuperCI>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
#if 0
  Dipole d(ref_->geom(), rdm1);
  d.compute();
#endif
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}
