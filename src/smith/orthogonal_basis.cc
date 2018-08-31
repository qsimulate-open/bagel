//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: orthogonal_basis.cc
// Copyright (C) 2018 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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
#ifdef COMPILE_SMITH

#include <numeric>
#include <src/smith/moint.h>
#include <src/smith/orthogonal_basis.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template<typename DataType>
Orthogonal_Basis<DataType>::Orthogonal_Basis(shared_ptr<const SMITH_Info<DataType>> info, const IndexRange c, const IndexRange a, const IndexRange v,
vector<double> f, vector<double> e0, shared_ptr<const Matrix> fact, shared_ptr<const Denom<DataType>> d, const bool residual,
shared_ptr<Vec<Tensor_<DataType>>> g0, shared_ptr<Vec<Tensor_<DataType>>> g1, shared_ptr<Vec<Tensor_<DataType>>> g2,
shared_ptr<Vec<Tensor_<DataType>>> g3, shared_ptr<Vec<Tensor_<DataType>>> g4) : closed_(c), active_(a), virt_(v), eig_(f), e0all_(e0), fockact_(fact) {
  nact_ = info->nact();
  nclosed_ = info->nclosed();
  nvirt_ = info->nvirt();
  nocc_ = nact_ + nclosed_;
  ncore_ = info->ncore();
  nclo_ = nclosed_ - ncore_;
  norb_ = nocc_ + nvirt_;
  nstates_ = info->ciwfn()->nstates();

  sssr_ = info->sssr();

  set_size(d);

  const int max = info->maxtile();

  interm_.push_back(IndexRange(0));
  for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
    interm_.push_back(IndexRange(shalf[iext]->ndim(), max));
  }

  for (int istate = 0; istate != nstates_; ++istate) {
    for (int ist = 0; ist != nstates_; ++ist) {
      auto tmp = make_shared<MultiTensor_<DataType>>(nstates_);
      if (!sssr_ || istate == ist) {
        for (int iext = Excitations::aibj; iext != Excitations::total; ++iext) {
          (*tmp)[iext + ist * Excitations::total] = init_data(iext);
        }
      }
      data_.push_back(tmp);
      denom_.push_back(tmp->clone());
    }
  }

  // let me store denominator for all orthogonal functions
  // TODO should we merge denom with it??
  set_denom(d);

  rdm0all_ = g0;
  rdm1all_ = g1;
  rdm2all_ = g2;
  rdm3all_ = g3;
  rdm4all_ = g4;

  if (residual) basis_type_ = Basis_Type::residual;
  else basis_type_ = Basis_Type::amplitude;

  zero();
}


template<typename DataType>
Orthogonal_Basis<DataType>::Orthogonal_Basis(const Orthogonal_Basis<DataType>& o, const bool clone, const bool residual) {
  closed_ = o.closed_;
  active_ = o.active_;
  virt_ = o.virt_;
  eig_ = o.eig_;
  e0all_ = o.e0all_;
  fockact_ = o.fockact_;

  nact_ = o.nact_;
  nclosed_ = o.nclosed_;
  nvirt_ = o.nvirt_;
  nocc_ = o.nocc_;
  ncore_ = o.ncore_;
  nclo_ = o.nclo_;
  norb_ = o.norb_;
  nstates_ = o.nstates_;

  sssr_ = o.sssr_;

  shalf_ = o.shalf_;
  size_ = o.size_;
  interm_ = o.interm_;

  data_ = o.data_;
  denom_ = o.denom_;

  rdm0all_ = o.rdm0all_;
  rdm1all_ = o.rdm1all_;
  rdm2all_ = o.rdm2all_;
  rdm3all_ = o.rdm3all_;
  rdm4all_ = o.rdm4all_;

  if (residual) basis_type_ = Basis_Type::residual;
  else basis_type_ = Basis_Type::amplitude;

  if (clone) zero();
}


template<typename DataType>
void Orthogonal_Basis<DataType>::set_size(shared_ptr<const Denom<DataType>> d) {
  const size_t size_aibj = nvirt_ * nvirt_ * nclo_ * nclo_;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt_ * nvirt_;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt_ * nclo_ * nvirt_;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo_ * nvirt_ * nclo_;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo_ * nclo_;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo_ * nvirt_;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt_;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo_;
  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;

  shalf_.push_back(make_shared<MatType>());
  shalf_.push_back(d->shalf_xx());
  shalf_.push_back(d->shalf_x());
  shalf_.push_back(d->shalf_h());
  shalf_.push_back(d->shalf_hh());
  shalf_.push_back(d->shalf_xh());
  shalf_.push_back(d->shalf_xxh());
  shalf_.push_back(d->shalf_xhh());

  size_.push_back(size_aibj);
  size_.push_back(size_arbs);
  size_.push_back(size_arbi);
  size_.push_back(size_airj);
  size_.push_back(size_risj);
  size_.push_back(size_airs);
  size_.push_back(size_arst);
  size_.push_back(size_rist);
  size_.push_back(size_all);
}


template<typename DataType>
shared_ptr<Tensor_<DataType>> Orthogonal_Basis<DataType>::init_data(const int iext) {
  // Now we use intermediate indices instead of orbital, and is somewhat complicated...
  // Interm always runs faster.
  unordered_set<size_t> sparse;
  shared_ptr<Tensor_<DataType>> out;
  switch(iext) {
    case Excitations::aibj:
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              sparse.insert(generate_hash_key(i0, i1, i2, i3));
      out = make_shared<Tensor_<DataType>>({closed_, virt_, closed_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::arbs:
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          for (auto& i0o : interm_[iext])
            sparse.insert(generate_hash_key(i0o, i1, i3));
      out = make_shared<Tensor_<DataType>>({interm_[iext], virt_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::arbi:
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext])
              sparse.insert(generate_hash_key(i0o, i1, i2, i3));
      out = make_shared<Tensor_<DataType>>({interm_[iext], virt_, closed_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::airj:
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext])
              sparse.insert(generate_hash_key(i0o, i0, i1, i2));
      out = make_shared<Tensor_<DataType>>({interm_[iext], closed_, virt_, closed_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::risj:
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext])
            sparse.insert(generate_hash_key(i0o, i0, i2));
      out = make_shared<Tensor_<DataType>>({interm_[iext], closed_, closed_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::airs:
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext])
            sparse.insert(generate_hash_key(i0o, i0, i1));
      out = make_shared<Tensor_<DataType>>({interm_[iext], closed_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::arst:
      for (auto& i1 : virt_)
        for (auto& i0o : interm_[iext])
          sparse.insert(generate_hash_key(i0o, i1));
      out = make_shared<Tensor_<DataType>>({interm_[iext], virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::rist:
      for (auto& i0 : closed_)
        for (auto& i0o : interm_[iext])
          sparse.insert(generate_hash_key(i0o, i0));
      out = make_shared<Tensor_<DataType>>({interm_[iext], closed_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
  }
  return out;
}


template<typename DataType>
void Orthogonal_Basis<DataType>::set_denom(shared_ptr<const Denom<DataType>> d) {
  for (int istate = 0; istate != nstates_; ++istate) {
    for (int ist = 0; ist != nstates_; ++ist) {
      e0_ = e0all_[istate];
      if (!sssr_ || ist == istate) {
        const double e0loc = e0all_[ist] - e0_;
        for (int iext = Excitations::aibj; iext != Excitations::total; ++iext) {
          const size_t interm_size = (iext == Excitations::aibj) ? 0 : shalf_[iext]->ndim();
          const shared_ptr<Tensor_<DataType>> dtensor = denom_[istate]->at(iext + ist * Excitations::total);
          switch(iext) {
            case Excitations::aibj:
              for (auto& i3 : virt_)
                for (auto& i2 : closed_)
                  for (auto& i1 : virt_)
                    for (auto& i0 : closed_) {
                      if (!dtensor->is_local(i0, i1, i2, i3)) continue;
                      unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0, i1, i2, i3)]);
                      size_t iall = 0;
                      for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                        for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                              data0[iall] = - eig_[j0] - eig_[j2] + eig_[j1] + eig_[j3] + e0loc;
                      denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0, i1, i2, i3);
                    }
              break;
            case Excitations::arbs:
              for (auto& i3 : virt_)
                for (auto& i1 : virt_)
                  for (auto& i0o : interm_[iext]) {
                    if (!dtensor->is_local(i0o, i1, i3)) continue;
                    unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i1, i3)]);
                    size_t iall = 0;
                    for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                      for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                        for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                          data0[iall] = eig_[j3] + eig_[j1] + d->denom_xx[j0o] - e0_;
                    denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i1, i3);
                  }
              break;
            case Excitations::arbi:
              for (auto& i3 : virt_)
                for (auto& i2 : closed_)
                  for (auto& i1 : virt_)
                    for (auto& i0o : interm_[iext]) {
                      if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
                      unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i1, i2, i3)]);
                      size_t iall = 0;
                      for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                        for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                            for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                              data0[iall] = eig_[j3] - eig_[j2] + eig_[j1] + d->denom_x[j0o] - e0_;
                      denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i1, i2, i3);
                    }
              break;
            case Excitations::airj:
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_)
                    for (auto& i0o : interm_[iext]) {
                      if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                      unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i0, i1, i2)]);
                      size_t iall = 0;
                      for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                        for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                          for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                            for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                              data0[iall] = - eig_[j2] + eig_[j1] - eig_[j0] + d->denom_h[j0o] - e0_;
                      denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i0, i1, i2);
                    }
              break;
            case Excitations::risj:
              for (auto& i2 : closed_)
                for (auto& i0 : closed_)
                  for (auto& i0o : interm_[iext]) {
                    if (!dtensor->is_local(i0o, i0, i2)) continue;
                    unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i0, i2)]);
                    size_t iall = 0;
                    for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                        for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                          data0[iall] = - eig_[j2] - eig_[j0] + d->denom_hh[j0o] - e0_;
                    denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i0, i2);
                  }
              break;
            case Excitations::airs:
              for (auto& i1 : virt_)
                for (auto& i0 : closed_)
                  for (auto& i0o : interm_[iext]) {
                    if (!dtensor->is_local(i0o, i0, i1)) continue;
                    unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i0, i1)]);
                    size_t iall = 0;
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                        for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                          data0[iall] = eig_[j1] - eig_[j0] + d->denom_xh[j0o] - e0_;
                    denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i0, i1);
                  }
              break;
            case Excitations::arst:
              for (auto& i1 : virt_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i1)) continue;
                  unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i1)]);
                  size_t iall = 0;
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                      data0[iall] = eig_[j1] + d->denom_xxh[j0o] - e0_;
                  denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i1);
                }
              break;
            case Excitations::rist:
              for (auto& i0 : closed_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i0)) continue;
                  unique_ptr<DataType[]> data0(new DataType[dtensor->get_size(i0o, i0)]);
                  size_t iall = 0;
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                      data0[iall] = - eig_[j0] + d->denom_xhh[j0o] - e0_;
                  denom_[istate]->at(iext + ist * Excitations::total)->put_block(data0, i0o, i0);
                }
              break;
          }
        }
      }
    }
  }
}


template<typename DataType>
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>,shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> Orthogonal_Basis<DataType>::feed_rdm(const int ist, const int jst) const {
  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  shared_ptr<RDM<3>> rdm3;
  shared_ptr<RDM<4>> rdm4;

  // collect den1ci
  {
    vector<IndexRange> o = rdm1all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    auto d1 = make_shared<RDM<1>>(nact_);
    for (auto& i1 : o[1].range())
      for (auto& i0 : o[0].range()) {
        auto input = rdm1all_->at(jst, ist)->get_block(i0, i1);
        for (size_t io1 = 0; io1 != i1.size(); ++io1)
          copy_n(&input[0 + i0.size() * io1], i0.size(), d1->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1));
      }
    rdm1 = d1->copy();
  }

  // collect den2ci
  {
    vector<IndexRange>o = rdm2all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    auto d2 = make_shared<RDM<2>>(nact_);
    for (auto& i3 : o[3].range())
      for (auto& i2 : o[2].range())
        for (auto& i1 : o[1].range())
          for (auto& i0 : o[0].range()) {
            auto input = rdm2all_->at(jst, ist)->get_block(i0, i1, i2, i3);
            for (size_t io3 = 0; io3 != i3.size(); ++io3)
              for (size_t io2 = 0; io2 != i2.size(); ++io2)
                for (size_t io1 = 0; io1 != i1.size(); ++io1)
                  copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * io3))], i0.size(),
                         d2->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2, io3 + i3.offset() - off3));
          }
    rdm2 = d2->copy();
  }

  // collect den3ci
  {
    vector<IndexRange>o = rdm3all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    const int off4 = o[4].front().offset();
    const int off5 = o[5].front().offset();
    auto d3 = make_shared<RDM<3>>(nact_);
    for (auto& i5 : o[5].range())
      for (auto& i4 : o[4].range())
        for (auto& i3 : o[3].range())
          for (auto& i2 : o[2].range())
            for (auto& i1 : o[1].range())
              for (auto& i0 : o[0].range()) {
                auto input = rdm3all_->at(jst, ist)->get_block(i0, i1, i2, i3, i4, i5);
                for (size_t io5 = 0; io5 != i5.size(); ++io5)
                  for (size_t io4 = 0; io4 != i4.size(); ++io4)
                    for (size_t io3 = 0; io3 != i3.size(); ++io3)
                      for (size_t io2 = 0; io2 != i2.size(); ++io2)
                        for (size_t io1 = 0; io1 != i1.size(); ++io1)
                          copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * io5))))],
                                 i0.size(), d3->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                 io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5));
              }
    rdm3 = d3->copy();
  }

  // collect den4ci
  {
    vector<IndexRange>o = rdm4all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    const int off4 = o[4].front().offset();
    const int off5 = o[5].front().offset();
    const int off6 = o[6].front().offset();
    const int off7 = o[7].front().offset();
    auto d4 = make_shared<RDM<4>>(nact_);
    for (auto& i7 : o[7].range())
      for (auto& i6 : o[6].range())
        for (auto& i5 : o[5].range())
          for (auto& i4 : o[4].range())
            for (auto& i3 : o[3].range())
              for (auto& i2 : o[2].range())
                for (auto& i1 : o[1].range())
                  for (auto& i0 : o[0].range()) {
                    auto input = rdm4all_->at(jst, ist)->get_block(i0, i1, i2, i3, i4, i5, i6, i7);
                    for (size_t io7 = 0; io7 != i7.size(); ++io7)
                      for (size_t io6 = 0; io6 != i6.size(); ++io6)
                        for (size_t io5 = 0; io5 != i5.size(); ++io5)
                          for (size_t io4 = 0; io4 != i4.size(); ++io4)
                            for (size_t io3 = 0; io3 != i3.size(); ++io3)
                              for (size_t io2 = 0; io2 != i2.size(); ++io2)
                                for (size_t io1 = 0; io1 != i1.size(); ++io1)
                                  copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * (io5 + i5.size() * (io6 + i6.size() * io7))))))],
                                         i0.size(), d4->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                         io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5, io6 + i6.offset() - off6, io7 + i7.offset() - off7));
              }
    rdm4 = d4->copy();
  }

  return tie(rdm1, rdm2, rdm3, rdm4);
}


template<typename DataType>
void Orthogonal_Basis<DataType>::transform_to_orthogonal(shared_ptr<const MultiTensor_<DataType>> t, const int istate) {
  // we put the transformed data in data_[istate].
  for (int ist = 0; ist != nstates_; ++ist) {
    e0_ = e0all_[istate];
    if (!sssr_ || ist == istate) {
      for (int iext = Excitations::aibj; iext != Excitations::total; ++iext) {
        const shared_ptr<Tensor_<DataType>> tensor = t->at(ist);
        const shared_ptr<Tensor_<DataType>> dtensor = denom_[istate]->at(iext + ist * Excitations::total);
        if (!dtensor) continue;
        switch(iext) {
          case Excitations::aibj:
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!tensor->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<DataType[]> data0 = tensor->get_block(i0, i1, i2, i3);
                    data_[istate]->at(iext + ist * Excitations::total)->add_block(data0, i0, i1, i2, i3);
                  }
            break;
          case Excitations::arbs:
            for (auto& i2 : active_)
              for (auto& i0 : active_) {
                auto create_transp = [this](const int i, const Index& I0, const Index& I2, const Index& I0o) {
                  unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*I0o.size()]);
                  for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
                    for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                      copy_n(denom_->shalf_xx()->element_ptr(I0o.offset(),(j0-nclosed_)+(j2-nclosed_)*nact_ + i*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*k);
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<DataType[]> transp = create_transp(ist, i0, i2, i0o);
                  for (auto& i3 : virt_)
                    for (auto& i1 : virt_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i0, i1, i2, i3);
                      unique_ptr<DataType[]> data1(new DataType[blocksize]);
                      // i0 i2 to be contracted
                      sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
                      unique_ptr<DataType[]> interm(new DataType[i1.size()*i3.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i1.size()*i3.size(), i0.size()*i2.size(),
                                                  sqrt(0.5), transp.get(), i0o.size(), data1.get(), i0.size()*i2.size(), 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i1, i3);
                    }
                }
              }
            break;
          case Excitations::arbi:
            for (auto& i0 : active_) {
              auto create_transp = [this](const int i, const Index& I0, const Index& I0o) {
                unique_ptr<DataType[]> out(new DataType[I0.size()*I0o.size()]);
                for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
                  copy_n(denom_->shalf_x()->element_ptr(I0o.offset(), j0-nclosed_ + i*nact_), I0o.size(), out.get()+I0o.size()*k);
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<DataType[]> transp = create_transp(ist, i0, i0o);
                for (auto& i3 : virt_)
                  for (auto& i2 : closed_)
                    for (auto& i1 : virt_) {
                      if (!tensor->is_local(i2, i3, i0, i1)) continue;
                      const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data2(new DataType[blocksize]);
                      // i0 to be contracted
                      sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
                      unique_ptr<DataType[]> interm(new DataType[i1.size()*i2.size()*i3.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i1.size()*i2.size()*i3.size(), i0.size(),
                                                  1.0, transp.get(), i0o.size(), data2.get(), i0.size(), 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i1, i2, i3);
                    }
              }
            }
            break;
          case Excitations::airj:
            for (auto& i3 : active_) {
              auto create_transp = [this](const int i, const Index& I3, const Index& I0o) {
                unique_ptr<DataType[]> out(new DataType[I3.size()*I0o.size()]);
                for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
                  copy_n(denom_->shalf_h()->element_ptr(I0o.offset(),j3-nclosed_ + i*nact_), I0o.size(), out.get()+I0o.size()*k);
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<DataType[]> transp = create_transp(ist, i3, i0o);
                for (auto& i2 : closed_)
                  for (auto& i1 : virt_)
                    for (auto& i0 : closed_) {
                      if (!tensor->is_local(i2, i3, i0, i1)) continue;
                      const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data2(new DataType[blocksize]);
                      // i3 to be contracted
                      sort_indices<1,2,3,0,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
                      unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*i2.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i0.size()*i1.size()*i2.size(), i3.size(),
                                                  1.0, transp.get(), i0o.size(), data2.get(), i3.size(), 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i0, i1, i2);
                    }
              }
            }
            break;
          case Excitations::risj:
            for (auto& i3 : active_)
              for (auto& i1 : active_) {
                auto create_transp = [this](const int i, const Index& I1, const Index& I3, const Index& I0o) {
                  unique_ptr<DataType[]> out(new DataType[I1.size()*I3.size()*I0o.size()]);
                  for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                    for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
                      copy_n(denom_->shalf_hh()->element_ptr(I0o.offset(),(j1-nclosed_)+(j3-nclosed_)*nact_ + i*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*k);
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<DataType[]> transp = create_transp(ist, i1, i3, i0o);
                  for (auto& i2 : closed_)
                    for (auto& i0 : closed_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i0, i1, i2, i3);
                      unique_ptr<DataType[]> data1(new DataType[blocksize]);
                      // i1 i3 to be contracted
                      sort_indices<1,3,0,2,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
                      unique_ptr<DataType[]> interm(new DataType[i0.size()*i2.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i0.size()*i2.size(), i1.size()*i3.size(),
                                                  sqrt(0.5), transp.get(), i0o.size(), data1.get(), i1.size()*i3.size(), 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i0, i2);
                    }
                }
              }
            break;
          case Excitations::airs:
            for (auto& i3 : active_)
              for (auto& i2 : active_) {
                auto create_transp = [this](const int i, const Index& I2, const Index& I3, const Index& I0o) {
                  unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*I0o.size()*2]);
                  for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                    for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
                      copy_n(denom_->shalf_xh()->element_ptr(I0o.offset(), (j2-nclosed_)+(j3-nclosed_)*nact_ + 2*i*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*k);
                      copy_n(denom_->shalf_xh()->element_ptr(I0o.offset(), (j2-nclosed_)+(j3-nclosed_)*nact_ + (2*i+1)*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*(k+I2.size()*I3.size()));
                    }
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<DataType[]> transp = create_transp(ist, i2, i3, i0o);
                  for (auto& i1 : virt_)
                    for (auto& i0 : closed_) {
                      if (!tensor->is_local(i2, i3, i0, i1)) continue;
                      const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data1 = tensor->get_block(i0, i3, i2, i1);
                      unique_ptr<DataType[]> data2(new DataType[blocksize*2]);
                      // i2 i3 to be contracted
                      sort_indices<0,1,2,3,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
                      sort_indices<2,1,0,3,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
                      unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i0.size()*i1.size(), i2.size()*i3.size()*2,
                                                  1.0, transp.get(), i0o.size(), data2.get(), i2.size()*i3.size()*2, 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i0, i1);
                    }
                }
              }
            break;
          case Excitations::arst:
            for (auto& i3 : active_)
              for (auto& i2 : active_)
                for (auto& i0 : active_) {
                  auto create_transp = [this](const int i, const Index& I0, const Index& I2, const Index& I3, const Index& I0o) {
                    unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*I3.size()*I0o.size()]);
                    for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                      for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
                        for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                          copy_n(denom_->shalf_xxh()->element_ptr(I0o.offset(),j0-nclosed_+nact_*(j2-nclosed_+nact_*(j3-nclosed_)) + i*nact_*nact_*nact_),
                                 I0o.size(), out.get()+I0o.size()*k);
                    return move(out);
                  };
                  for (auto& i0o : interm_[iext]) {
                    unique_ptr<DataType[]> transp = create_transp(ist, i0, i2, i3, i0o);
                    for (auto& i1 : virt_) {
                      if (!tensor->is_local(i2, i3, i0, i1)) continue;
                      const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data1(new DataType[blocksize]);
                      // i3 i2 i0 to be contracted
                      sort_indices<1,0,2,3,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
                      unique_ptr<DataType[]> interm(new DataType[i1.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i1.size(), i0.size()*i2.size()*i3.size(),
                                                  1.0, transp.get(), i0o.size(), data1.get(), i0.size()*i2.size()*i3.size(), 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i1);
                    }
                  }
                }
            break;
          case Excitations::rist:
            for (auto& i3 : active_)
              for (auto& i1 : active_)
                for (auto& i0 : active_) {
                  auto create_transp = [this](const int i, const Index& I0, const Index& I1, const Index& I3, const Index& I0o) {
                    unique_ptr<DataType[]> out(new DataType[I0.size()*I1.size()*I3.size()*I0o.size()]);
                    for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                      for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
                        for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                          copy_n(denom_->shalf_xhh()->element_ptr(I0o.offset(),j0-nclosed_+nact_*(j1-nclosed_+nact_*(j3-nclosed_)) + i*nact_*nact_*nact_),
                                 I0o.size(), out.get()+I0o.size()*k);
                    return move(out);
                  };
                  for (auto& i0o : interm_[iext]) {
                    unique_ptr<DataType[]> transp = create_transp(ist, i0, i1, i3, i0o);
                    for (auto& i2 : closed_) {
                      if (!tensor->is_local(i2, i3, i0, i1)) continue;
                      const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
                      unique_ptr<DataType[]> data1(new DataType[blocksize]);
                      // i3 i1 i0 to be contracted
                      sort_indices<1,3,2,0,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
                      unique_ptr<DataType[]> interm(new DataType[i2.size()*i0o.size()]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i2.size(), i0.size()*i1.size()*i3.size(),
                                                  1.0, transp.get(), i0o.size(), data1.get(), i0.size()*i1.size()*i3.size(), 0.0, interm.get(), i0o.size());
                      denom_[istate]->at(iext + ist * Excitations::total)->add_block(interm, i0o, i2);
                    }
                  }
                }
            break;
        }
      }
    }
  }
}


template<typename DataType>
vector<shared_ptr<MultiTensor_<DataType>>> Orthogonal_Basis<DataType>::transform_to_redundant() {
  vector<shared_ptr<MultiTensor_<DataType>>> out;

  for (int istate = 0; istate != nstates_; ++istate)
    out.push_back(transform_to_redundant(istate));

  return out;
}


template<typename DataType>
shared_ptr<MultiTensor_<DataType>> Orthogonal_Basis<DataType>::transform_to_redundant(const int istate) {
  auto out = make_shared<MultiTensor_<DataType>>();
  return out;
}


template<typename DataType>
void Orthogonal_Basis<DataType>::update(shared_ptr<const Orthogonal_Basis<DataType>> r, const double shift, const bool imag) {
}


template<typename DataType>
void Orthogonal_Basis<DataType>::print_convergence(shared_ptr<const Orthogonal_Basis<DataType>> s, shared_ptr<const Orthogonal_Basis<DataType>> r) {
}

#endif
