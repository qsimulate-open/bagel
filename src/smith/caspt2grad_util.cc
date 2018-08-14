//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: caspt2grad_util.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <bagel_config.h>
#include <src/scf/hf/fock.h>
#include <src/smith/caspt2grad.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Matrix> CASPT2Grad::diagonal_D1() const {
#ifdef COMPILE_SMITH
  vector<IndexRange> ind = d2_->indexrange();
  if (!(ind[0] == ind[2] && ind[1] == ind[3])) throw logic_error("wrong");

  Tensor interm(vector<IndexRange>{ind[0], ind[3]});
  interm.allocate();

  for (auto& i3 : ind[3])
    for (auto& i1 : ind[1])
      for (auto& i0 : ind[0]) {
        if (!interm.is_local(i0, i3) || !d2_->exists(i0, i1, i1, i3)) continue;
        unique_ptr<double[]> in = d2_->get_block(i0, i1, i1, i3);
        unique_ptr<double[]> buf(new double[i0.size()*i3.size()]);
        fill_n(buf.get(), i0.size()*i3.size(), 0.0);
        for (int j3 = 0; j3 != i3.size(); ++j3)
          for (int j1 = 0; j1 != i1.size(); ++j1)
            blas::ax_plus_y_n(1.0, in.get()+i0.size()*(j1+i1.size()*(j1+i1.size()*j3)),
                              i0.size(), buf.get()+i0.size()*j3);
        interm.add_block(buf, i0, i3);
      }
  return interm.matrix();
#else
  return coeff_->clone(); // dummy
#endif
}


shared_ptr<Matrix> CASPT2Grad::spin_density_unrelaxed() const {
#ifdef COMPILE_SMITH
  const int nele_act = fci_->det()->nelea() + fci_->det()->neleb();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nmo = d11_->mdim();

  const double sp1 = fci_->det()->nspin()*0.5 + 1.0;

  auto out1 = make_shared<Matrix>(*d11_);
  auto out0 = out1->clone();
  shared_ptr<const Matrix> d0 = ref_->rdm1_mat(target_);
  out0->add_block(1.0, 0, 0, nocc, nocc, d0);
  out0->scale((4.0 - nele_act) * 0.5);
  out1->scale((4.0 - nele_act) * 0.5);

  // Correcting for first order RDMs
  // two-body part
  out1->add_block(-2.0, ncore_, nclosed, nocc-ncore_, nmo-nclosed, diagonal_D1());
  // closed part
  out1->add_block(-4.0, ncore_, nclosed, nclosed-ncore_, nmo-nclosed, d11_->get_submatrix(ncore_, nclosed, nclosed-ncore_, nmo-nclosed));

  // Correcting for zeroth order RDMs
  // closed part
  out0->add_block(nele_act*0.5-2.0, 0, 0, nclosed, nclosed, d0->get_submatrix(0, 0, nclosed, nclosed));
  // two-body part
  shared_ptr<const RDM<2>> rdm2 = ref_->rdm2(target());
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        (*out0)(j+nclosed,i+nclosed) -= rdm2->element(j,k,k,i);

  // scale with 1/(S+1)
  out0->scale(1.0 / sp1);
  out1->scale(1.0 / sp1);
  // finally symmetrize
  out0->symmetrize();
  out1->symmetrize();

  // add second-order contribution
  assert(target() == 0); // this is assumed in CASPT2 when computing the norm
  // CAUTION - this is because spin-free d(2) includes the <1|1>d(0) term, but the alpha spin one does not.
  shared_ptr<const Matrix> d0active = ref_->rdm1(target_)->rdm1_mat(nclosed, false)->resize(nmo,nmo);
  auto out2 = make_shared<Matrix>(*sd1_ - (*out0 + *d0active) * wf1norm_[target()]);

  if (fabs(out2->trace()) > 1.0e-8)
    cout << "  **** warning **** Trace of the second-order CASPT2 spin-density matrix is nonzero: "
         << setprecision(10) << out2->trace() << endl;
  if ((*out1 - *sd11_).rms() > 1.0e-8)
    cout << "  **** warning **** First order CASPT2 spin-density matrix has inconsistency: "
         << (*out1 - *sd11_).rms() << endl;

  // returns an AO matrix
  return make_shared<Matrix>(*coeff_ * (*out0 + *out1 + *out2) ^ *coeff_);
#else
  return coeff_->clone(); // dummy
#endif
}


// relaxation part of the spin density
shared_ptr<Matrix> CASPT2Grad::spin_density_relax(shared_ptr<const RDM<1>> zrdm1, shared_ptr<const RDM<2>> zrdm2, shared_ptr<const Matrix> zmat) const {
#ifdef COMPILE_SMITH
  // zrdm1 and 2 are defined only withtin the active space
  const int nele_act = fci_->det()->nelea() + fci_->det()->neleb();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nmo  = nocc + ref_->nvirt();

  shared_ptr<const RDM<1>> rdm1 = ref_->rdm1_av();
  shared_ptr<const RDM<2>> rdm2 = ref_->rdm2_av();
  Matrix rdmsa(nact, nact);

  auto out = make_shared<Matrix>(nmo, nmo);

  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j) {
      (*out)(j+nclosed, i+nclosed) = (4 - nele_act) * 0.5 * zrdm1->element(j, i);
      rdmsa(j, i) = (4 - nele_act) * 0.5 * rdm1->element(j, i);
      for (int k = 0; k != nact; ++k) {
        (*out)(j+nclosed, i+nclosed) -= zrdm2->element(j, k, k, i);
        rdmsa(j, i) -= rdm2->element(j, k, k, i);
      }
    }
  // add to the block
  out->add_block(2.0, 0, nclosed, nmo, nact, (zmat->slice(nclosed, nocc) * rdmsa));

  // scale with 1/(S+1) that is common in both contributions
  out->scale(1.0 / (fci_->det()->nspin()*0.5 + 1.0));
  // finally symmetrize
  out->symmetrize();
  return make_shared<Matrix>(*coeff_ * *out ^ *coeff_);
#else
  return zmat->clone(); // dummy
#endif
}


shared_ptr<DFFullDist> CASPT2Grad::contract_D1(shared_ptr<const DFFullDist> full) const {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nocc = ref_->nocc();
  const int nall = nocc + ref_->nvirt();

  shared_ptr<DFFullDist> out = full->clone();
  double* optr = out->block(0)->data();
  const double* iptr = full->block(0)->data();
  const size_t asize = full->block(0)->asize();

  // first create a tensor for full.
  // Aux index blocking
  shared_ptr<const StaticDist> adist = ref_->geom()->df()->adist_now();
  const IndexRange aux(adist);

  vector<IndexRange> ind = d2_->indexrange();
  assert(ind[0] == ind[2] && ind[1] == ind[3]);

  // create a GA array
  auto input = make_shared<Tensor>(vector<IndexRange>{aux, ind[0], ind[1]});
  input->allocate();
  auto output = input->clone();

  // fill into input
  for (auto& a : aux) {
    tuple<size_t, size_t> loc = adist->locate(a.offset());
    if (get<0>(loc) != mpi__->rank())
      continue;
    for (auto& j : ind[1])
      for (auto& i : ind[0]) {
        unique_ptr<double[]> buf(new double[a.size()*i.size()*j.size()]);
        for (int ej = 0; ej != j.size(); ++ej)
          for (int ei = 0; ei != i.size(); ++ei)
            copy_n(iptr + get<1>(loc)+asize*(i.offset()+ei+full->nocc1()*(j.offset()+ej)), a.size(), buf.get()+a.size()*(ei+i.size()*ej));
        input->put_block(buf, a, i, j);
      }
  }
  mpi__->barrier();

  // next contract
  for (auto& j : ind[1])
    for (auto& i : ind[0])
      for (auto& a : aux) {
        if (!output->is_local(a, i, j)) continue;
        unique_ptr<double[]> buf(new double[a.size()*i.size()*j.size()]);
        fill_n(buf.get(), a.size()*i.size()*j.size(), 0.0);
        for (auto& l : ind[1]) {
          for (auto& k : ind[0]) {
            unique_ptr<double[]> in = input->get_block(a, k, l);
            if (d2_->exists(k, l, i, j)) {
              unique_ptr<double[]> d = d2_->get_block(k, l, i, j);
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, a.size(), i.size()*j.size(), k.size()*l.size(),
                                          1.0, in.get(), a.size(), d.get(), k.size()*l.size(), 1.0, buf.get(), a.size());
            } else if (d2_->exists(i, j, k, l)) {
              unique_ptr<double[]> d = d2_->get_block(i, j, k, l);
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, a.size(), i.size()*j.size(), k.size()*l.size(),
                                          1.0, in.get(), a.size(), d.get(), i.size()*j.size(), 1.0, buf.get(), a.size());
            }
          }
        }
        output->add_block(buf, a, i, j);
      }
  mpi__->barrier();

  // finally fill back into output
  for (auto& a : aux) {
    tuple<size_t, size_t> loc = adist->locate(a.offset());
    if (get<0>(loc) != mpi__->rank())
      continue;
    for (auto& j : ind[1])
      for (auto& i : ind[0]) {
        unique_ptr<double[]> buf = output->get_block(a, i, j);
        for (int ej = 0; ej != j.size(); ++ej)
          for (int ei = 0; ei != i.size(); ++ei)
            copy_n(buf.get()+a.size()*(ei+i.size()*ej), a.size(), optr + get<1>(loc)+asize*(i.offset()+ei+full->nocc1()*(j.offset()+ej)));
      }
  }
  mpi__->barrier();

  // one body part
  {
    auto is_cc   = [&](const int& i) { return i < nclosed; };
    auto is_act  = [&](const int& i) { return i >= nclosed && i < nocc; };
    auto is_virt = [&](const int& i) { return i >= nocc; };
    auto is_exc  = [&](const int& i, const int& j) { return is_virt(i) || (is_act(i) && is_cc(j)); };

    for (int t = 0; t != nall; ++t)
      for (int l = 0; l != nocc; ++l)
        // c(cc)x, c(cc)a, x(cc)a
        if (is_exc(t, l))
          for (int s = 0; s != nclosed; ++s) {
            blas::ax_plus_y_n(2.0*d11_->element(l,t), iptr + asize*(l+nocc*t), asize, optr + asize*(s+nocc*s));
            blas::ax_plus_y_n(   -d11_->element(l,t), iptr + asize*(s+nocc*t), asize, optr + asize*(l+nocc*s));
            blas::ax_plus_y_n(2.0*d11_->element(l,t), iptr + asize*(s+nocc*s), asize, optr + asize*(l+nocc*t));
            blas::ax_plus_y_n(   -d11_->element(l,t), iptr + asize*(l+nocc*s), asize, optr + asize*(s+nocc*t));
          }
  }
  return out;
#else
  return full->clone(); // dummy
#endif
}


void CASPT2Grad::augment_Y(shared_ptr<Matrix> d0ms, shared_ptr<Matrix> g0, shared_ptr<Dvec> g1, shared_ptr<const DFHalfDist> halfj, const int istate, const int jstate, const double egap) {
#ifdef COMPILE_SMITH
  const int nmobasis = coeff_->mdim();

  g0->add_block(egap, 0, 0, nmobasis, nmobasis, *vd1_);

  for (int ist = 0; ist != nstates_; ++ist) {
    for (int jst = 0; jst != nstates_; ++jst) {
      if (ist != jst) {
        const double fac = .5 * (msrot(ist, istate) * msrot(jst, jstate) - msrot(ist, jstate) * msrot(jst, istate)) * egap;
        g1->data(jst)->ax_plus_y(fac, fci_->civectors()->data(ist));
      }
    }
  }
#endif
}


tuple<shared_ptr<Matrix>, shared_ptr<const DFFullDist>>
  CASPT2Grad::compute_Y(shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> halfj, shared_ptr<const DFHalfDist> halfjj, const bool nacme) {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nmobasis = coeff_->mdim();

  const MatView ocmat = coeff_->slice(0, nocc);

  shared_ptr<const DFFullDist> full = halfj->compute_second_transform(coeff_);
  shared_ptr<const DFFullDist> fullo = halfj->compute_second_transform(ocmat);

  // Y_rs = 2[Y1 + Y2 + Y3(ri) + Y4 + Y5(ri)]
  shared_ptr<Matrix> out = make_shared<Matrix>(nmobasis, nmobasis);
  const MatView ocoeff = coeff_->slice(0, nocc);

  shared_ptr<RDM<1>> rdm10;
  shared_ptr<RDM<2>> rdm20;
  {
    // 2 Y1 = h(d0 + d1 + d2) * 2
    // one-electron contributions
    const Matrix hmo(*coeff_ % *ref_->hcore() * *coeff_);
    auto d0 = make_shared<Matrix>(nmobasis, nmobasis);
    if (nact)
      d0->add_block(1.0, nclosed, nclosed, nact, nact, *d10ms_);
    if (nacme) {
      d0->symmetrize();
    } else {
      for (int i = 0; i != nclosed; ++i)
        d0->element(i,i) = 2.0;
    }
    *out += hmo * (*d1_ + *d11_ + *d0) * 2.0;
  }

  {
    // Y2 = Y2_rs = Y2_ri + Y2_ra, so making both at once
    shared_ptr<Matrix> dkl;
    dkl = ref_->rdm1_mat(); // D0SA
    dkl->sqrt();
    dkl->scale(1.0/sqrt(2.0));
    Fock<1> fock(geom_, ref_->hcore()->clone(), nullptr, ocoeff * *dkl, /*grad*/false, /*rhf*/true);
    *out += *coeff_ % fock * *coeff_ * *d1_ * 2.0;
  }

  {
    // 2 Y3 = 2 Y3_ri*dm0_ji
    // coulomb
    auto dmrao = make_shared<Matrix>(*coeff_ * *d1_ ^ *coeff_);
    auto jop = geom_->df()->compute_Jop(dmrao);
    // exchange
    auto kopi = halfjj->compute_Kop_1occ(dmrao, -0.5)->transpose();

    auto tmp = make_shared<Matrix>(*coeff_ % (*jop * ocoeff + *kopi));
    if (nact) {
      *tmp *= *ref_->rdm1_mat(); // D0SA
    } else {
      *tmp *= 2.0;
    }
    out->add_block(2.0, 0, 0, nmobasis, nocc, tmp);
  }

  // fullks will be reused for gradient contraction
  shared_ptr<const DFFullDist> fullks;
  {
    // 2 Y4 =  2 K^{kl}_{rt} D^{lk}_{ts} = 2 (kr|lj) D0_(lj,ki) +  2 (kr|lt) D1_(lt,ks)
    // construct stepwise, D1 part
    fullks = contract_D1(full);
    *out += *full->form_2index(fullks, 2.0);
    // D0 part
    if (nacme) {
      shared_ptr<const DFFullDist> fulld = fullo->apply_2rdm_tran(*d20ms_, *d10ms_, nclosed, nact);
      out->add_block(2.0, 0, 0, nmobasis, nocc, full->form_2index(fulld, 0.5));
      out->add_block(2.0, 0, 0, nmobasis, nocc, full->form_2index(fulld->swap(), 0.5));
    } else {
      shared_ptr<const DFFullDist> fulld = nact ? fullo->apply_2rdm(*d20ms_, *d10ms_, nclosed, nact) : fullo->apply_closed_2RDM();
      out->add_block(2.0, 0, 0, nmobasis, nocc, full->form_2index(fulld, 1.0));
    }
  }

  {
    // 2 Y5 = 2 Y5_ri = 2 Ybar (Gyorffy)  = 2 (rs|tj) D^ij_st = 2 (rl|jk) D0_(il,jk) + 2 (rs|tj) D1_(is,jt)]
    // construct stepwise, D1 part
    shared_ptr<const DFHalfDist> dfback = fullks->apply_J()->back_transform(coeff_);
    auto y5ri_ao = ref_->geom()->df()->form_2index(dfback, 1.0);
    out->add_block(2.0, 0, 0, nmobasis, nocc, *coeff_ % *y5ri_ao);
  }

  return make_tuple(out, fullks);
#else
  return tuple<shared_ptr<Matrix>, shared_ptr<const DFFullDist>>();
#endif
}
