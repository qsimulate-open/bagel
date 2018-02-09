//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci/distfci.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <src/ci/fci/distfci.h>

using namespace std;
using namespace bagel;


void DistFCI::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_->at(0).
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norb_);
    rdm2_av_ = make_shared<RDM<2>>(norb_);
  } else if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  // we need expanded lists
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, /*compressed=*/false, /*mute=*/true);
  cc_->set_det(detex);

  for (int i = 0; i != nstate_; ++i)
    compute_rdm12(i, i);

  // calculate state averaged RDMs
  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weight_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weight_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0,0);
    rdm2_av_ = rdm2_->at(0,0);
  }

  cc_->set_det(det_);
}


void DistFCI::compute_rdm12(const int ist, const int jst) {
  if (det_->compress()) {
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }

  shared_ptr<DistCivec> ccbra = cc_->data(ist);
  shared_ptr<DistCivec> ccket = cc_->data(jst);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_civec(ccbra, ccket);

  // setting to private members.
  rdm1_->emplace(ist, jst, rdm1);
  rdm2_->emplace(ist, jst, rdm2);

  cc_->set_det(det_);
}

tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> DistFCI::rdm34(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>();
}


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<3>>> DistFCI::rdm34f(const int ist, const int jst, shared_ptr<const Matrix> fock) const {
  return tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<3>>>();
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> DistFCI::rdm12_alpha(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>();
}


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> DistFCI::rdm34_alpha(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>();
}


// calculate <I|a+a|0>
void DistFCI::sigma_2a1(shared_ptr<const DistCivec> cc, shared_ptr<DistDvec> d) const {
  assert(d->det() == cc->det());
  const int lb = cc->lenb();
  for (int ip = 0; ip != d->ij(); ++ip) {
    shared_ptr<DistCivec> tcc = d->data(ip);
    for (auto& i : cc->det()->phia(ip))
      if (tcc->is_local(i.source)) {
        unique_ptr<double[]> source = cc->rma_get(i.target);
        const double sign = static_cast<double>(i.sign);
        blas::scale_n(sign, source.get(), lb);
        tcc->rma_add(source, i.source);
      }
  }
}

// calculate <I|b+b|0>
void DistFCI::sigma_2a2(shared_ptr<const DistCivec> cc, shared_ptr<DistDvec> d) const {
  assert(d->det() == cc->det());
  const int lb = cc->lenb();
  for (int i = 0; i != cc->lena(); ++i) {
    if (!cc->is_local(i)) continue;
    unique_ptr<double[]> source = cc->rma_get(i);
    unique_ptr<double[]> target(new double[lb]);
    for (int ip = 0; ip != d->ij(); ++ip) {
      fill_n(target.get(), lb, 0.0);
      for (auto& iter : cc->det()->phib(ip))
        target[iter.source] += iter.sign * source[iter.target];
      assert(d->data(ip)->is_local(i));
      d->data(ip)->rma_add(target, i);
    }
  }
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
DistFCI::compute_rdm12_from_civec(shared_ptr<const DistCivec> cbra, shared_ptr<const DistCivec> cket) const {

  // since we consider here number conserving operators...
  auto dbra = make_shared<DistDvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<DistDvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<DistDvec>(cket->det(), norb_*norb_);
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }
  return compute_rdm12_last_step(dbra, dket, cbra);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
DistFCI::compute_rdm12_last_step(shared_ptr<const DistDvec> dbra, shared_ptr<const DistDvec> dket, shared_ptr<const DistCivec> cibra) const {

  const int nri = cibra->asize()*cibra->lenb();
  const int ij  = norb_*norb_;

  // 1RDM c^dagger <I|\hat{E}|0>
  // 2RDM \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, ij);
    for (int i = 0; i != ij; ++i)
      copy_n(dket->data(i)->data(), nri, dket_data->element_ptr(0, i));
    auto rdm1t = btas::group(*rdm1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, rdm1t, {1});

    auto dbra_data = dket_data;
    if (dbra != dket) {
      dbra_data = make_shared<Matrix>(nri, ij);
      for (int i = 0; i != ij; ++i)
        copy_n(dbra->data(i)->data(), nri, dbra_data->element_ptr(0, i));
    }
    auto rdm2t = group(group(*rdm2, 2,4), 0,2);
    btas::contract(1.0, *dbra_data, {1,0}, *dket_data, {1,2}, 0.0, rdm2t, {0,2});
  }

  rdm1->allreduce();
  rdm2->allreduce();

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, rdm2->element_ptr(0,0,k,i));
    }
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb_; ++i)
    for (int k = 0; k != norb_; ++k)
      for (int j = 0; j != norb_; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}


shared_ptr<Dvec> DistFCI::rdm1deriv(const int istate) const {

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<DistCivec> cbra = cc_->data(istate);

  // 1RDM ci derivative
  // <I|E_ij|0>

  auto dbra = make_shared<DistDvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  return distdvec_to_dvec(dbra);
}


shared_ptr<Dvec> DistFCI::rdm2deriv(const int istate) const {
  throw logic_error("DistFCI::rdm2deriv is not implemented yet...");
  return nullptr;
}


shared_ptr<Matrix> DistFCI::rdm2fderiv(const int istate, shared_ptr<const Matrix> fock, shared_ptr<const Matrix> dmat) const {
  throw logic_error("DistFCI::rdm2fderiv is not implemented yet...");
  return nullptr;
}


shared_ptr<Matrix> DistFCI::rdm2deriv_offset(const int istate, const size_t dsize, const size_t offset, shared_ptr<const Matrix> dmat, const bool parallel) const {
  throw logic_error("DistFCI::rdm2deriv_offset is not implemented yet...");
  return nullptr;
}


tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>
DistFCI::rdm3deriv(const int istate, shared_ptr<const Matrix> fock, const size_t offset, const size_t size, shared_ptr<const Matrix> dbra_in, shared_ptr<const Matrix> fock_ebra_in) const {
  throw logic_error("DistFCI::rdm3deriv is not implemented yet...");
  return tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>();
}

