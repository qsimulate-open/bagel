//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: orthogonal.cc
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
#include <src/smith/orthogonal.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

Orthogonal_Basis::Orthogonal_Basis(shared_ptr<const SMITH_Info<double>> info, const IndexRange c, const IndexRange a, const IndexRange v,
vector<double> f, vector<double> e0, shared_ptr<const Matrix> fact, shared_ptr<const Denom<double>> d, const bool residual,
shared_ptr<Vec<Tensor_<double>>> g0, shared_ptr<Vec<Tensor_<double>>> g1, shared_ptr<Vec<Tensor_<double>>> g2,
shared_ptr<Vec<Tensor_<double>>> g3, shared_ptr<Vec<Tensor_<double>>> g4) : closed_(c), active_(a), virt_(v), eig_(f), e0all_(e0) {
  nact_ = info->nact();
  nclosed_ = info->nclosed();
  nvirt_ = info->nvirt();
  nocc_ = nact_ + nclosed_;
  ncore_ = info->ncore();
  nclo_ = nclosed_ - ncore_;
  norb_ = nocc_ + nvirt_;
  nstates_ = info->ciwfn()->nstates();

  fockact_ = fact->copy();
  sssr_ = info->sssr();
  imag_ = info->shift_imag();
  shift_ = info->shift();

  set_size(d);

  const int max = info->maxtile();

  int keyoffset = closed_.nblock() + active_.nblock() + virt_.nblock();
  for (int iext = Excitations::arbs; iext != Excitations::aibj; ++iext) {
    interm_.push_back(IndexRange(shalf_[iext]->ndim(), max, keyoffset));
    keyoffset += interm_[iext].nblock();
  }
  interm_.push_back(IndexRange(0));

//  cout << " istate" << endl;
  for (int istate = 0; istate != nstates_; ++istate) {
    auto tmp = make_shared<MultiTensor_<double>>(Excitations::total + (sssr_ ? 0 : nstates_-1));
    for (int iext = Excitations::arbs; iext != Excitations::aibj; ++iext) {
      (*tmp)[iext] = init_data(iext);
    }
    // aibj depends on whether sssr or msmr
    if (sssr_) {
      const int pos = Excitations::aibj;
      (*tmp)[pos] = init_data(Excitations::aibj);
    } else {
      for (int ist = 0; ist != nstates_; ++ist) {
        const int pos = Excitations::aibj + ist;
        (*tmp)[pos] = init_data(Excitations::aibj);
      }
    }

    data_.push_back(tmp);
    denom_.push_back(tmp->clone());
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


Orthogonal_Basis::Orthogonal_Basis(const Orthogonal_Basis& o, const bool clone, const bool residual) {
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
  imag_ = o.imag_;
  shift_ = o.shift_;

  shalf_.resize(Excitations::total);
  for (int i = 0; i != Excitations::total; ++i)
    shalf_[i] = o.shalf_[i];

  size_ = o.size_;
  interm_ = o.interm_;

  data_.resize(nstates_);
  for (int i = 0; i != nstates_; ++i) {
    if (clone) {
      data_[i] = o.data_[i]->clone();
    } else {
      data_[i] = o.data_[i]->copy();
    }
  }

  denom_.resize(nstates_);
  for (int i = 0; i != nstates_; ++i)
    denom_[i] = o.denom_[i]->copy();

  rdm0all_ = o.rdm0all_;
  rdm1all_ = o.rdm1all_;
  rdm2all_ = o.rdm2all_;
  rdm3all_ = o.rdm3all_;
  rdm4all_ = o.rdm4all_;

  if (residual) basis_type_ = Basis_Type::residual;
  else basis_type_ = Basis_Type::amplitude;
}


void Orthogonal_Basis::set_size(shared_ptr<const Denom<double>> d) {
  const size_t size_arbs = d->shalf_xx()->ndim()  * nvirt_ * nvirt_;
  const size_t size_arbi = d->shalf_x()->ndim()   * nvirt_ * nclo_ * nvirt_;
  const size_t size_airj = d->shalf_h()->ndim()   * nclo_ * nvirt_ * nclo_;
  const size_t size_risj = d->shalf_hh()->ndim()  * nclo_ * nclo_;
  const size_t size_airs = d->shalf_xh()->ndim()  * nclo_ * nvirt_;
  const size_t size_arst = d->shalf_xxh()->ndim() * nvirt_;
  const size_t size_rist = d->shalf_xhh()->ndim() * nclo_;
  const size_t size_aibj = nvirt_ * nvirt_ * nclo_ * nclo_ * (sssr_ ? 1 : nstates_);
  const size_t size_all = size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist + size_aibj;

  shalf_.push_back(d->shalf_xx()->copy());
  shalf_.push_back(d->shalf_x()->copy());
  shalf_.push_back(d->shalf_h()->copy());
  shalf_.push_back(d->shalf_hh()->copy());
  shalf_.push_back(d->shalf_xh()->copy());
  shalf_.push_back(d->shalf_xxh()->copy());
  shalf_.push_back(d->shalf_xhh()->copy());
  shalf_.push_back(make_shared<Matrix>());

  size_.push_back(size_arbs);
  size_.push_back(size_arbi);
  size_.push_back(size_airj);
  size_.push_back(size_risj);
  size_.push_back(size_airs);
  size_.push_back(size_arst);
  size_.push_back(size_rist);
  size_.push_back(size_aibj);
  size_.push_back(size_all);
}


shared_ptr<Tensor_<double>> Orthogonal_Basis::init_data(const int iext) {
  // Now we use intermediate indices instead of orbital, and is somewhat complicated...
  // Interm always runs faster.
  unordered_set<size_t> sparse;
  shared_ptr<Tensor_<double>> out;
  switch(iext) {
    case Excitations::arbs:
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          for (auto& i0o : interm_[iext])
            sparse.insert(generate_hash_key(i0o, i1, i3));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], virt_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::arbi:
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext])
              sparse.insert(generate_hash_key(i0o, i1, i2, i3));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], virt_, closed_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::airj:
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext])
              sparse.insert(generate_hash_key(i0o, i0, i1, i2));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], closed_, virt_, closed_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::risj:
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext])
            sparse.insert(generate_hash_key(i0o, i0, i2));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], closed_, closed_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::airs:
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext])
            sparse.insert(generate_hash_key(i0o, i0, i1));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], closed_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::arst:
      for (auto& i1 : virt_)
        for (auto& i0o : interm_[iext])
          sparse.insert(generate_hash_key(i0o, i1));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::rist:
      for (auto& i0 : closed_)
        for (auto& i0o : interm_[iext])
          sparse.insert(generate_hash_key(i0o, i0));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{interm_[iext], closed_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
    case Excitations::aibj:
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              sparse.insert(generate_hash_key(i0, i1, i2, i3));
      out = make_shared<Tensor_<double>>(vector<IndexRange>{closed_, virt_, closed_, virt_}, /*kramers=*/false, sparse, /*alloc=*/true);
      break;
  }
  return out;
}


shared_ptr<MultiTensor_<double>> Orthogonal_Basis::weight_by_denom(const int istate, shared_ptr<const MultiTensor_<double>> original) const {
  auto out = make_shared<MultiTensor_<double>>(Excitations::total + (sssr_ ? 0 : nstates_-1));
  for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
    const shared_ptr<Tensor_<double>> dtensor = denom_[istate]->at(iext);
    shared_ptr<const Tensor_<double>> ttensor = original->at(iext);
    switch(iext) {
      case Excitations::arbs:
        out->at(iext) = ttensor->clone();
        for (auto& i3 : virt_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i1, i3)) continue;
              const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i1, i3);
              unique_ptr<double[]> tdata = ttensor->get_block(i0o, i1, i3);
              const size_t intermsize = dtensor->get_size(i0o, i1, i3);
              for (size_t i = 0; i != intermsize; ++i) {
                tdata[i] *= ddata[i];
              }
              out->at(iext)->add_block(tdata, i0o, i1, i3);
            }
        break;
      case Excitations::arbi:
        out->at(iext) = ttensor->clone();
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
                const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i1, i2, i3);
                unique_ptr<double[]> tdata = ttensor->get_block(i0o, i1, i2, i3);
                const size_t intermsize = dtensor->get_size(i0o, i1, i2, i3);
                for (size_t i = 0; i != intermsize; ++i) {
                  tdata[i] *= ddata[i];
                }
                out->at(iext)->add_block(tdata, i0o, i1, i2, i3);
              }
        break;
      case Excitations::airj:
        out->at(iext) = ttensor->clone();
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i0, i1, i2);
                unique_ptr<double[]> tdata = ttensor->get_block(i0o, i0, i1, i2);
                const size_t intermsize = dtensor->get_size(i0o, i0, i1, i2);
                for (size_t i = 0; i != intermsize; ++i) {
                  tdata[i] *= ddata[i];
                }
                out->at(iext)->add_block(tdata, i0o, i0, i1, i2);
              }
        break;
      case Excitations::risj:
        out->at(iext) = ttensor->clone();
        for (auto& i2 : closed_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i0, i2)) continue;
              const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i0, i2);
              unique_ptr<double[]> tdata = ttensor->get_block(i0o, i0, i2);
              const size_t intermsize = dtensor->get_size(i0o, i0, i2);
              for (size_t i = 0; i != intermsize; ++i) {
                tdata[i] *= ddata[i];
              }
              out->at(iext)->add_block(tdata, i0o, i0, i2);
            }
        break;
      case Excitations::airs:
        out->at(iext) = ttensor->clone();
        for (auto& i1 : virt_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i0, i1)) continue;
              const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i0, i1);
              unique_ptr<double[]> tdata = ttensor->get_block(i0o, i0, i1);
              const size_t intermsize = dtensor->get_size(i0o, i0, i1);
              for (size_t i = 0; i != intermsize; ++i) {
                tdata[i] *= ddata[i];
              }
              out->at(iext)->add_block(tdata, i0o, i0, i1);
            }
        break;
      case Excitations::arst:
        out->at(iext) = ttensor->clone();
        for (auto& i1 : virt_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i1)) continue;
            const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i1);
            unique_ptr<double[]> tdata = ttensor->get_block(i0o, i1);
            const size_t intermsize = dtensor->get_size(i0o, i1);
            for (size_t i = 0; i != intermsize; ++i) {
              tdata[i] *= ddata[i];
            }
            out->at(iext)->add_block(tdata, i0o, i1);
          }
        break;
      case Excitations::rist:
        out->at(iext) = ttensor->clone();
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i0)) continue;
            const unique_ptr<double[]> ddata = dtensor->get_block(i0o, i0);
            unique_ptr<double[]> tdata = ttensor->get_block(i0o, i0);
            const size_t intermsize = dtensor->get_size(i0o, i0);
            for (size_t i = 0; i != intermsize; ++i) {
              tdata[i] *= ddata[i];
            }
            out->at(iext)->add_block(tdata, i0o, i0);
          }
       break;
      case Excitations::aibj:
        for (int ist = 0; ist != nstates_; ++ist) {
          if (!sssr_ || ist == istate) {
            const int pos = iext + (sssr_ ? 0 : ist);
            out->at(pos) = original->at(pos)->clone();
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!denom_[istate]->at(pos)->is_local(i0, i1, i2, i3)) continue;
                      const unique_ptr<double[]> ddata = denom_[istate]->at(pos)->get_block(i0, i1, i2, i3);
                      unique_ptr<double[]> tdata = original->at(pos)->get_block(i0, i1, i2, i3);
                      const size_t intermsize = denom_[istate]->at(pos)->get_size(i0, i1, i2, i3);
                      for (size_t i = 0; i != intermsize; ++i) {
                        tdata[i] *= ddata[i];
                      }
                      out->at(pos)->add_block(tdata, i0, i1, i2, i3);
                  }
          }
        }
      break;
    }
  }
  mpi__->barrier();
  return out;
}


shared_ptr<MultiTensor_<double>> Orthogonal_Basis::get_contravariant(const int istate, const bool weight) const {
  auto out = make_shared<MultiTensor_<double>>(Excitations::total + (sssr_ ? 0 : nstates_-1));
  for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
    const shared_ptr<Tensor_<double>> dtensor = data_[istate]->at(iext);
    switch(iext) {
      case Excitations::arbs:
        out->at(iext) = dtensor->copy();
        break;
      case Excitations::arbi:
        out->at(iext) = dtensor->clone();
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
                unique_ptr<double[]> data0 = dtensor->get_block(i0o, i1, i2, i3);
                unique_ptr<double[]> data1 = dtensor->get_block(i0o, i3, i2, i1);
                sort_indices<0,3,2,1,2,1,-1,1>(data1, data0, i0o.size(), i3.size(), i2.size(), i1.size());
                out->at(iext)->add_block(data0, i0o, i1, i2, i3);
              }
        break;
      case Excitations::airj:
        out->at(iext) = dtensor->clone();
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                unique_ptr<double[]> data0 = dtensor->get_block(i0o, i0, i1, i2);
                unique_ptr<double[]> data1 = dtensor->get_block(i0o, i2, i1, i0);
                sort_indices<0,3,2,1,2,1,-1,1>(data1, data0, i0o.size(), i2.size(), i1.size(), i0.size());
                out->at(iext)->add_block(data0, i0o, i0, i1, i2);
              }
        break;
      case Excitations::risj:
        out->at(iext) = dtensor->copy();
        break;
      case Excitations::airs:
        out->at(iext) = dtensor->copy();
        break;
      case Excitations::arst:
        out->at(iext) = dtensor->copy();
        break;
      case Excitations::rist:
        out->at(iext) = dtensor->copy();
        break;
      case Excitations::aibj:
        for (int ist = 0; ist != nstates_; ++ist) {
          if (!sssr_ || ist == istate) {
            const int pos = iext + (sssr_ ? 0 : ist);
            out->at(pos) = data_[istate]->at(pos)->clone();
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!data_[istate]->at(pos)->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<double[]> data0 = data_[istate]->at(pos)->get_block(i0, i1, i2, i3);
                    unique_ptr<double[]> data1 = data_[istate]->at(pos)->get_block(i0, i3, i2, i1);
                    sort_indices<0,3,2,1,8,1,-4,1>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
                    out->at(pos)->add_block(data0, i0, i1, i2, i3);
                  }
          }
        }
        break;
    }
  }
  mpi__->barrier();

  if (weight && imag_) {
    // weight by denominator, so that we can get the correct energy
    out = weight_by_denom(istate, out);
  }
  return out;
}


void Orthogonal_Basis::set_denom(shared_ptr<const Denom<double>> d) {
  for (int istate = 0; istate != nstates_; ++istate) {
    e0_ = e0all_[istate];
    for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
      const shared_ptr<Tensor_<double>> dtensor = denom_[istate]->at(iext);
      switch(iext) {
        case Excitations::arbs:
          for (auto& i3 : virt_)
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i1, i3)) continue;
                unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i1, i3)]);
                size_t iall = 0;
                for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                      data0[iall] = eig_[j3] + eig_[j1] + d->denom_xx(j0o) - e0_;
                denom_[istate]->at(iext)->put_block(data0, i0o, i1, i3);
              }
          break;
        case Excitations::arbi:
          for (auto& i3 : virt_)
            for (auto& i2 : closed_)
              for (auto& i1 : virt_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
                  unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i1, i2, i3)]);
                  size_t iall = 0;
                  for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                    for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                      for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                        for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                          data0[iall] = eig_[j3] - eig_[j2] + eig_[j1] + d->denom_x(j0o) - e0_;
                  denom_[istate]->at(iext)->put_block(data0, i0o, i1, i2, i3);
                }
          break;
        case Excitations::airj:
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0 : closed_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                  unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i0, i1, i2)]);
                  size_t iall = 0;
                  for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                        for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                          data0[iall] = - eig_[j2] + eig_[j1] - eig_[j0] + d->denom_h(j0o) - e0_;
                  denom_[istate]->at(iext)->put_block(data0, i0o, i0, i1, i2);
                }
          break;
        case Excitations::risj:
          for (auto& i2 : closed_)
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i0, i2)) continue;
                unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i0, i2)]);
                size_t iall = 0;
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                      data0[iall] = - eig_[j2] - eig_[j0] + d->denom_hh(j0o) - e0_;
                denom_[istate]->at(iext)->put_block(data0, i0o, i0, i2);
              }
          break;
        case Excitations::airs:
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i0, i1)) continue;
                unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i0, i1)]);
                size_t iall = 0;
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                      data0[iall] = eig_[j1] - eig_[j0] + d->denom_xh(j0o) - e0_;
                denom_[istate]->at(iext)->put_block(data0, i0o, i0, i1);
              }
          break;
        case Excitations::arst:
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i1)) continue;
              unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i1)]);
              size_t iall = 0;
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                  data0[iall] = eig_[j1] + d->denom_xxh(j0o) - e0_;
              denom_[istate]->at(iext)->put_block(data0, i0o, i1);
            }
          break;
        case Excitations::rist:
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i0)) continue;
              unique_ptr<double[]> data0(new double[dtensor->get_size(i0o, i0)]);
              size_t iall = 0;
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall)
                  data0[iall] = - eig_[j0] + d->denom_xhh(j0o) - e0_;
              denom_[istate]->at(iext)->put_block(data0, i0o, i0);
            }
         break;
        case Excitations::aibj:
          for (int ist = 0; ist != nstates_; ++ist) {
            double e0loc = e0_ - e0all_[ist];
            if (!sssr_ || ist == istate) {
              const int pos = iext + (sssr_ ? 0 : ist);
              for (auto& i3 : virt_)
                for (auto& i2 : closed_)
                  for (auto& i1 : virt_)
                    for (auto& i0 : closed_) {
                      if (!denom_[istate]->at(pos)->is_local(i0, i1, i2, i3)) continue;
                      unique_ptr<double[]> data0(new double[denom_[istate]->at(pos)->get_size(i0, i1, i2, i3)]);
                      size_t iall = 0;
                      for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                        for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                              data0[iall] = - eig_[j0] - eig_[j2] + eig_[j1] + eig_[j3] + e0loc;
                      denom_[istate]->at(pos)->put_block(data0, i0, i1, i2, i3);
                    }
              }
            }
            break;
      }
    }
  }
  mpi__->barrier();
}


void Orthogonal_Basis::transform_to_orthogonal(shared_ptr<const MultiTensor_<double>> t, int istate) {
  // we put the transformed data in data_[istate].
  for (int ist = 0; ist != nstates_; ++ist) {
    if (!t->at(ist)) continue;
    for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
      shared_ptr<const Tensor_<double>> tensor = t->at(ist);
      switch(iext) {
        case Excitations::arbs:
          for (auto& i2 : active_)
            for (auto& i0 : active_) {
              auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I2, const Index& I0o) {
                unique_ptr<double[]> out(new double[I0.size()*I2.size()*I0o.size()]);
                for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
                  for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                    copy_n(shalf_[iext]->element_ptr(I0o.offset(),(j0-nclosed_)+(j2-nclosed_)*nact_ + i*nact_*nact_),
                           I0o.size(), out.get()+I0o.size()*k);
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<double[]> transp = create_transp(ist, i0, i2, i0o);
                for (auto& i3 : virt_)
                  for (auto& i1 : virt_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i1, i3)) continue;
                    const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                    const unique_ptr<double[]> data0 = tensor->get_block(i0, i1, i2, i3);
                    unique_ptr<double[]> data1(new double[blocksize]);
                    // i0 i2 to be contracted
                    sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
                    unique_ptr<double[]> interm(new double[i1.size()*i3.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i1.size()*i3.size(), i0.size()*i2.size(),
                                                sqrt(0.5), transp.get(), i0o.size(), data1.get(), i0.size()*i2.size(), 0.0, interm.get(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm, i0o, i1, i3);
                  }
              }
            }
          break;
        case Excitations::arbi:
          for (auto& i0 : active_) {
            auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I0o) {
              unique_ptr<double[]> out(new double[I0.size()*I0o.size()]);
              for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
                copy_n(shalf_[iext]->element_ptr(I0o.offset(), j0-nclosed_ + i*nact_), I0o.size(), out.get()+I0o.size()*k);
              return move(out);
            };
            for (auto& i0o : interm_[iext]) {
              unique_ptr<double[]> transp = create_transp(ist, i0, i0o);
              for (auto& i3 : virt_)
                for (auto& i2 : closed_)
                  for (auto& i1 : virt_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i1, i2, i3)) continue;
                    const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                    const unique_ptr<double[]> data0 = tensor->get_block(i2, i3, i0, i1);
                    unique_ptr<double[]> data2(new double[blocksize]);
                    // i0 to be contracted
                    sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
                    const unique_ptr<double[]> data1 = tensor->get_block(i2, i1, i0, i3);
                    sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());
                    unique_ptr<double[]> interm(new double[i1.size()*i2.size()*i3.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i1.size()*i2.size()*i3.size(), i0.size(),
                                                1.0, transp.get(), i0o.size(), data2.get(), i0.size(), 0.0, interm.get(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm, i0o, i1, i2, i3);
                  }
            }
          }
          break;
        case Excitations::airj:
          for (auto& i3 : active_) {
            auto create_transp = [&iext, this](const int i, const Index& I3, const Index& I0o) {
              unique_ptr<double[]> out(new double[I3.size()*I0o.size()]);
              for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
                copy_n(shalf_[iext]->element_ptr(I0o.offset(),j3-nclosed_ + i*nact_), I0o.size(), out.get()+I0o.size()*k);
              return move(out);
            };
            for (auto& i0o : interm_[iext]) {
              unique_ptr<double[]> transp = create_transp(ist, i3, i0o);
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i0, i1, i2)) continue;
                    const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                    const unique_ptr<double[]> data0 = tensor->get_block(i2, i3, i0, i1);
                    unique_ptr<double[]> data2(new double[blocksize]);
                    // i3 to be contracted
                    sort_indices<1,2,3,0,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
                    const unique_ptr<double[]> data1 = tensor->get_block(i0, i3, i2, i1);
                    sort_indices<1,0,3,2,2,3,1,3>(data1, data2, i0.size(), i3.size(), i2.size(), i1.size());
                    unique_ptr<double[]> interm(new double[i0.size()*i1.size()*i2.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i0.size()*i1.size()*i2.size(), i3.size(),
                                                1.0, transp.get(), i0o.size(), data2.get(), i3.size(), 0.0, interm.get(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm, i0o, i0, i1, i2);
                  }
            }
          }
          break;
        case Excitations::risj:
          for (auto& i3 : active_)
            for (auto& i1 : active_) {
              auto create_transp = [&iext, this](const int i, const Index& I1, const Index& I3, const Index& I0o) {
                unique_ptr<double[]> out(new double[I1.size()*I3.size()*I0o.size()]);
                for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                  for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
                    copy_n(shalf_[iext]->element_ptr(I0o.offset(),(j1-nclosed_)+(j3-nclosed_)*nact_ + i*nact_*nact_),
                           I0o.size(), out.get()+I0o.size()*k);
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<double[]> transp = create_transp(ist, i1, i3, i0o);
                for (auto& i2 : closed_)
                  for (auto& i0 : closed_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i0, i2)) continue;
                    const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                    const unique_ptr<double[]> data0 = tensor->get_block(i0, i1, i2, i3);
                    unique_ptr<double[]> data1(new double[blocksize]);
                    // i1 i3 to be contracted
                    sort_indices<1,3,0,2,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
                    unique_ptr<double[]> interm(new double[i0.size()*i2.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i0.size()*i2.size(), i1.size()*i3.size(),
                                                sqrt(0.5), transp.get(), i0o.size(), data1.get(), i1.size()*i3.size(), 0.0, interm.get(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm, i0o, i0, i2);
                  }
              }
            }
          break;
        case Excitations::airs:
          for (auto& i3 : active_)
            for (auto& i2 : active_) {
              auto create_transp = [&iext, this](const int i, const Index& I2, const Index& I3, const Index& I0o) {
                unique_ptr<double[]> out(new double[I2.size()*I3.size()*I0o.size()*2]);
                for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                  for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
                    copy_n(shalf_[iext]->element_ptr(I0o.offset(), (j2-nclosed_)+(j3-nclosed_)*nact_ + 2*i*nact_*nact_),
                           I0o.size(), out.get()+I0o.size()*k);
                    copy_n(shalf_[iext]->element_ptr(I0o.offset(), (j2-nclosed_)+(j3-nclosed_)*nact_ + (2*i+1)*nact_*nact_),
                           I0o.size(), out.get()+I0o.size()*(k+I2.size()*I3.size()));
                  }
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<double[]> transp = create_transp(ist, i2, i3, i0o);
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i0, i1)) continue;
                    const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                    const unique_ptr<double[]> data0 = tensor->get_block(i2, i3, i0, i1);
                    const unique_ptr<double[]> data1 = tensor->get_block(i0, i3, i2, i1);
                    unique_ptr<double[]> data2(new double[blocksize*2]);
                    // i2 i3 to be contracted
                    sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
                    sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
                    unique_ptr<double[]> interm(new double[i0.size()*i1.size()*i0o.size()]);
                    unique_ptr<double[]> interm2(new double[i0.size()*i1.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size(), i0o.size(), i2.size()*i3.size()*2,
                                                1.0, data2.get(), i0.size()*i1.size(), transp.get(), i0o.size(), 0.0, interm.get(), i0.size()*i1.size());
                    sort_indices<2,0,1,0,1,1,1>(interm.get(), interm2.get(), i0.size(), i1.size(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm2, i0o, i0, i1);
                  }
              }
            }
          break;
        case Excitations::arst:
          for (auto& i3 : active_)
            for (auto& i2 : active_)
              for (auto& i0 : active_) {
                auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I2, const Index& I3, const Index& I0o) {
                  unique_ptr<double[]> out(new double[I0.size()*I2.size()*I3.size()*I0o.size()]);
                  for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                    for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
                      for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                        copy_n(shalf_[iext]->element_ptr(I0o.offset(),j0-nclosed_+nact_*(j2-nclosed_+nact_*(j3-nclosed_)) + i*nact_*nact_*nact_),
                               I0o.size(), out.get()+I0o.size()*k);
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<double[]> transp = create_transp(ist, i0, i2, i3, i0o);
                  for (auto& i1 : virt_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i1)) continue;
                    const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                    const unique_ptr<double[]> data0 = tensor->get_block(i2, i3, i0, i1);
                    unique_ptr<double[]> data1(new double[blocksize]);
                    // i3 i2 i0 to be contracted
                    sort_indices<2,0,1,3,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
                    unique_ptr<double[]> interm(new double[i1.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i1.size(), i0.size()*i2.size()*i3.size(),
                                                1.0, transp.get(), i0o.size(), data1.get(), i0.size()*i2.size()*i3.size(), 0.0, interm.get(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm, i0o, i1);
                  }
                }
              }
          break;
        case Excitations::rist:
          for (auto& i3 : active_)
            for (auto& i1 : active_)
              for (auto& i0 : active_) {
                auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I1, const Index& I3, const Index& I0o) {
                  unique_ptr<double[]> out(new double[I0.size()*I1.size()*I3.size()*I0o.size()]);
                  for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                    for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
                      for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                        copy_n(shalf_[iext]->element_ptr(I0o.offset(),j0-nclosed_+nact_*(j1-nclosed_+nact_*(j3-nclosed_)) + i*nact_*nact_*nact_),
                               I0o.size(), out.get()+I0o.size()*k);
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<double[]> transp = create_transp(ist, i0, i1, i3, i0o);
                  for (auto& i2 : closed_) {
                    if (!data_[istate]->at(iext)->is_local(i0o, i2)) continue;
                    const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
                    const unique_ptr<double[]> data0 = tensor->get_block(i2, i3, i0, i1);
                    unique_ptr<double[]> data1(new double[blocksize]);
                    // i3 i1 i0 to be contracted
                    sort_indices<2,3,1,0,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
                    unique_ptr<double[]> interm(new double[i2.size()*i0o.size()]);
                    btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0o.size(), i2.size(), i0.size()*i1.size()*i3.size(),
                                                1.0, transp.get(), i0o.size(), data1.get(), i0.size()*i1.size()*i3.size(), 0.0, interm.get(), i0o.size());
                    data_[istate]->at(iext)->add_block(interm, i0o, i2);
                  }
                }
              }
          break;
        case Excitations::aibj:
          if (!sssr_ || ist == istate) {
            const int pos = iext + (sssr_ ? 0 : ist);
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!data_[istate]->at(iext)->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<double[]> data0 = t->at(ist)->get_block(i0, i1, i2, i3);
                    const unique_ptr<double[]> data1 = t->at(ist)->get_block(i0, i3, i2, i1);
                    sort_indices<0,3,2,1,2,12,1,12>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
                    data_[istate]->at(pos)->add_block(data0, i0, i1, i2, i3);
                  }
          }
          break;
      }
    }
  }
  mpi__->barrier();
}


shared_ptr<MultiTensor_<double>> Orthogonal_Basis::transform_to_redundant(const int istate) const {
  // we put the transformed data in out.
  auto out = make_shared<MultiTensor_<double>>(nstates_);

  for (int ist = 0; ist != nstates_; ++ist) {
    if (!sssr_ || ist == istate) {
      (*out)[ist] = init_amplitude();
      for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
        shared_ptr<Tensor_<double>> tensor = out->at(ist);
        shared_ptr<const Tensor_<double>> dtensor = data_[istate]->at(iext);
        switch(iext) {
          case Excitations::arbs:
            for (auto& i2 : active_)
              for (auto& i0 : active_) {
                auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I2, const Index& I0o) {
                  unique_ptr<double[]> out(new double[I0.size()*I2.size()*I0o.size()]);
                  for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
                    for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                      copy_n(shalf_[iext]->element_ptr(I0o.offset(),(j0-nclosed_)+(j2-nclosed_)*nact_ + i*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*k);
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<double[]> transp = create_transp(ist, i0, i2, i0o);
                  for (auto& i3 : virt_)
                    for (auto& i1 : virt_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i1, i3);
                      const size_t interm_size = dtensor->get_size(i0o, i1, i3);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o, i1, i3);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> data0(new double[blocksize]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), i0o.size(),
                                                  sqrt(0.5), transp.get(), i0o.size(), interm.get(), i0o.size(), 0.0, data0.get(), i0.size()*i2.size());
                      unique_ptr<double[]> data1(new double[blocksize]);
                      sort_indices<0,2,1,3,0,1,1,1>(data0.get(), data1.get(), i0.size(), i2.size(), i1.size(), i3.size());
                      tensor->add_block(data1, i0, i1, i2, i3);
                    }
                }
              }
            break;
          case Excitations::arbi:
            for (auto& i0 : active_) {
              auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I0o) {
                unique_ptr<double[]> out(new double[I0.size()*I0o.size()]);
                for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
                  copy_n(shalf_[iext]->element_ptr(I0o.offset(), j0-nclosed_ + i*nact_), I0o.size(), out.get()+I0o.size()*k);
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<double[]> transp = create_transp(ist, i0, i0o);
                for (auto& i3 : virt_)
                  for (auto& i2 : closed_)
                    for (auto& i1 : virt_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i1, i2, i3);
                      const size_t interm_size = dtensor->get_size(i0o, i1, i2, i3);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o, i1, i2, i3);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> data0(new double[blocksize]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size(), i1.size()*i2.size()*i3.size(), i0o.size(),
                                                  1.0, transp.get(), i0o.size(), interm.get(), i0o.size(), 0.0, data0.get(), i0.size());
                      tensor->add_block(data0, i0, i1, i2, i3);
                    }
              }
            }
            break;
          case Excitations::airj:
            for (auto& i3 : active_) {
              auto create_transp = [&iext, this](const int i, const Index& I3, const Index& I0o) {
                unique_ptr<double[]> out(new double[I3.size()*I0o.size()]);
                for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
                  copy_n(shalf_[iext]->element_ptr(I0o.offset(),j3-nclosed_ + i*nact_), I0o.size(), out.get()+I0o.size()*k);
                return move(out);
              };
              for (auto& i0o : interm_[iext]) {
                unique_ptr<double[]> transp = create_transp(ist, i3, i0o);
                for (auto& i2 : closed_)
                  for (auto& i1 : virt_)
                    for (auto& i0 : closed_) {
                      if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i0, i1, i2);
                      const size_t interm_size = dtensor->get_size(i0o, i0, i1, i2);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o,i0, i1, i2);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> interm2(new double[interm_size]);
                      sort_indices<1,2,3,0,0,1,1,1>(interm.get(), interm2.get(), i0o.size(), i0.size(), i1.size(), i2.size());
                      unique_ptr<double[]> data0(new double[blocksize]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size()*i2.size(), i3.size(), i0o.size(),
                                                  1.0, interm2.get(), i0.size()*i1.size()*i2.size(), transp.get(), i0o.size(), 0.0, data0.get(), i0.size()*i1.size()*i2.size());
                      tensor->add_block(data0, i0, i1, i2, i3);
                    }
              }
            }
            break;
          case Excitations::risj:
            for (auto& i3 : active_)
              for (auto& i1 : active_) {
                auto create_transp = [&iext, this](const int i, const Index& I1, const Index& I3, const Index& I0o) {
                  unique_ptr<double[]> out(new double[I1.size()*I3.size()*I0o.size()]);
                  for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                    for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
                      copy_n(shalf_[iext]->element_ptr(I0o.offset(),(j1-nclosed_)+(j3-nclosed_)*nact_ + i*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*k);
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<double[]> transp = create_transp(ist, i1, i3, i0o);
                  for (auto& i2 : closed_)
                    for (auto& i0 : closed_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i0, i2);
                      const size_t interm_size = dtensor->get_size(i0o, i0, i2);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o, i0, i2);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> interm2(new double[interm_size]);
                      sort_indices<1,2,0,0,1,1,1>(interm.get(), interm2.get(), i0o.size(), i0.size(), i2.size());
                      unique_ptr<double[]> data0(new double[blocksize]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), i0o.size(),
                                                  sqrt(0.5), interm2.get(), i0.size()*i2.size(), transp.get(), i0o.size(), 0.0, data0.get(), i0.size()*i2.size());
                      unique_ptr<double[]> data1(new double[blocksize]);
                      sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
                      tensor->add_block(data1, i0, i1, i2, i3);
                    }
                }
              }
            break;
          case Excitations::airs:
            for (auto& i3 : active_)
              for (auto& i2 : active_) {
                auto create_transp = [&iext, this](const int i, const Index& I2, const Index& I3, const Index& I0o) {
                  unique_ptr<double[]> out(new double[I2.size()*I3.size()*I0o.size()*2]);
                  for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                    for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
                      copy_n(shalf_[iext]->element_ptr(I0o.offset(), (j2-nclosed_)+(j3-nclosed_)*nact_ + 2*i*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*k);
                      copy_n(shalf_[iext]->element_ptr(I0o.offset(), (j2-nclosed_)+(j3-nclosed_)*nact_ + (2*i+1)*nact_*nact_),
                             I0o.size(), out.get()+I0o.size()*(k+I2.size()*I3.size()));
                    }
                  return move(out);
                };
                for (auto& i0o : interm_[iext]) {
                  unique_ptr<double[]> transp = create_transp(ist, i2, i3, i0o);
                  for (auto& i1 : virt_)
                    for (auto& i0 : closed_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i0, i1);
                      const size_t interm_size = dtensor->get_size(i0o, i0, i1);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o, i0, i1);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> interm2(new double[dtensor->get_size(i0o, i0, i1)]);
                      sort_indices<1,2,0,0,1,1,1>(interm.get(), interm2.get(), i0o.size(), i0.size(), i1.size());
                      unique_ptr<double[]> data0(new double[blocksize*2]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2.size()*i3.size()*2, i0o.size(),
                                                  1.0, interm2.get(), i0.size()*i1.size(), transp.get(), i0o.size(), 0.0, data0.get(), i0.size()*i1.size());
                      unique_ptr<double[]> data1(new double[blocksize]);
                      unique_ptr<double[]> data2(new double[blocksize]);
                      copy_n(data0.get(), blocksize, data1.get());
                      sort_indices<2,1,0,3,0,1,1,1>(data0.get()+blocksize, data2.get(), i0.size(), i1.size(), i2.size(), i3.size());
                      tensor->add_block(data1, i0, i1, i2, i3);
                      tensor->add_block(data2, i2, i1, i0, i3);
                    }
                }
              }
            break;
          case Excitations::arst:
            for (auto& i3 : active_)
              for (auto& i2 : active_)
                for (auto& i0 : active_) {
                  auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I2, const Index& I3, const Index& I0o) {
                    unique_ptr<double[]> out(new double[I0.size()*I2.size()*I3.size()*I0o.size()]);
                    for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                      for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
                        for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                          copy_n(shalf_[iext]->element_ptr(I0o.offset(),j0-nclosed_+nact_*(j2-nclosed_+nact_*(j3-nclosed_)) + i*nact_*nact_*nact_),
                                 I0o.size(), out.get()+I0o.size()*k);
                    return move(out);
                  };
                  for (auto& i0o : interm_[iext]) {
                    unique_ptr<double[]> transp = create_transp(ist, i0, i2, i3, i0o);
                    for (auto& i1 : virt_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i1);
                      const size_t interm_size = dtensor->get_size(i0o, i1);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o, i1);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> data0(new double[blocksize]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size()*i2.size()*i3.size(), i1.size(), i0o.size(),
                                                  1.0, transp.get(), i0o.size(), interm.get(), i0o.size(), 0.0, data0.get(), i0.size()*i2.size()*i3.size());
                      unique_ptr<double[]> data1(new double[blocksize]);
                      sort_indices<0,3,1,2,0,1,1,1>(data0.get(), data1.get(), i0.size(), i2.size(), i3.size(), i1.size());
                      tensor->add_block(data1, i0, i1, i2, i3);
                    }
                  }
                }
            break;
          case Excitations::rist:
            for (auto& i3 : active_)
              for (auto& i1 : active_)
                for (auto& i0 : active_) {
                  auto create_transp = [&iext, this](const int i, const Index& I0, const Index& I1, const Index& I3, const Index& I0o) {
                    unique_ptr<double[]> out(new double[I0.size()*I1.size()*I3.size()*I0o.size()]);
                    for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
                      for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
                        for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                          copy_n(shalf_[iext]->element_ptr(I0o.offset(),j0-nclosed_+nact_*(j1-nclosed_+nact_*(j3-nclosed_)) + i*nact_*nact_*nact_),
                                 I0o.size(), out.get()+I0o.size()*k);
                    return move(out);
                  };
                  for (auto& i0o : interm_[iext]) {
                    unique_ptr<double[]> transp = create_transp(ist, i0, i1, i3, i0o);
                    for (auto& i2 : closed_) {
                      if (!tensor->is_local(i0, i1, i2, i3)) continue;
                      const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
                      unique_ptr<double[]> interm = dtensor->get_block(i0o, i2);
                      const size_t interm_size = dtensor->get_size(i0o, i2);
                      if (imag_) {
                        const unique_ptr<double[]> denom = denom_[istate]->at(iext)->get_block(i0o, i2);
                        for (size_t i = 0; i != interm_size; ++i) {
                          interm[i] *= denom[i];
                        }
                      }
                      unique_ptr<double[]> data0(new double[blocksize]);
                      btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size()*i1.size()*i3.size(), i2.size(), i0o.size(),
                                                  1.0, transp.get(), i0o.size(), interm.get(), i0o.size(), 0.0, data0.get(), i0.size()*i1.size()*i3.size());
                      unique_ptr<double[]> data1(new double[blocksize]);
                      sort_indices<0,1,3,2,0,1,1,1>(data0.get(), data1.get(), i0.size(), i1.size(), i3.size(), i2.size());
                      tensor->add_block(data1, i0, i1, i2, i3);
                    }
                  }
                }
            break;
          case Excitations::aibj:
            const int pos = iext + (sssr_ ? 0 : ist);
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!tensor->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<double[]> data0 = data_[istate]->at(pos)->get_block(i0, i1, i2, i3);
                    if (imag_) {
                      const unique_ptr<double[]> denom = denom_[istate]->at(pos)->get_block(i0,i1, i2, i3);
                      const size_t blocksize = data_[istate]->at(pos)->get_size(i0, i1, i2, i3);
                      for (size_t i = 0; i != blocksize; ++i) {
                        data0[i] *= denom[i];
                      }
                    }
                    out->at(ist)->add_block(data0, i0, i1, i2, i3);
                  }
            break;
        }
      }
    }
  }
  mpi__->barrier();

  return out;
}


void Orthogonal_Basis::update(shared_ptr<const Orthogonal_Basis> r, const int istate) {
  const double shift2 = shift_ * shift_;

  for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
    const shared_ptr<Tensor_<double>> rtensor = r->data(istate)->at(iext);
    const shared_ptr<Tensor_<double>> dtensor = denom_[istate]->at(iext);
    switch(iext) {
      case Excitations::arbs:
        for (auto& i3 : virt_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i1, i3)) continue;
              unique_ptr<double[]> residual = rtensor->get_block(i0o, i1, i3);
              unique_ptr<double[]> denom    = dtensor->get_block(i0o, i1, i3);
              const size_t blocksize = rtensor->get_size(i0o, i1, i3);
              if (imag_) {
                for (size_t j = 0; j != blocksize; ++j) {
                  residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
                }
              } else {
                for (size_t j = 0; j != blocksize; ++j) {
                  residual[j] *= -(1.0 / (denom[j] + shift_));
                }
              }
              data_[istate]->at(iext)->add_block(residual, i0o, i1, i3);
            }
        break;
      case Excitations::arbi:
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
                unique_ptr<double[]> residual = rtensor->get_block(i0o, i1, i2, i3);
                unique_ptr<double[]> denom    = dtensor->get_block(i0o, i1, i2, i3);
                const size_t blocksize = rtensor->get_size(i0o, i1, i2, i3);
                if (imag_) {
                  for (size_t j = 0; j != blocksize; ++j) {
                    residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
                  }
                } else {
                  for (size_t j = 0; j != blocksize; ++j) {
                    residual[j] *= -(1.0 / (denom[j] + shift_));
                  }
                }
                data_[istate]->at(iext)->add_block(residual, i0o, i1, i2, i3);
              }
        break;
      case Excitations::airj:
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                unique_ptr<double[]> residual = rtensor->get_block(i0o, i0, i1, i2);
                unique_ptr<double[]> denom    = dtensor->get_block(i0o, i0, i1, i2);
                const size_t blocksize = rtensor->get_size(i0o, i0, i1, i2);
                if (imag_) {
                  for (size_t j = 0; j != blocksize; ++j) {
                    residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
                  }
                } else {
                  for (size_t j = 0; j != blocksize; ++j) {
                    residual[j] *= -(1.0 / (denom[j] + shift_));
                  }
                }
                data_[istate]->at(iext)->add_block(residual, i0o, i0, i1, i2);
              }
        break;
      case Excitations::risj:
        for (auto& i2 : closed_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i0, i2)) continue;
              unique_ptr<double[]> residual = rtensor->get_block(i0o, i0, i2);
              unique_ptr<double[]> denom    = dtensor->get_block(i0o, i0, i2);
              const size_t blocksize = rtensor->get_size(i0o, i0, i2);
              if (imag_) {
                for (size_t j = 0; j != blocksize; ++j) {
                  residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
                }
              } else {
                for (size_t j = 0; j != blocksize; ++j) {
                  residual[j] *= -(1.0 / (denom[j] + shift_));
                }
              }
              data_[istate]->at(iext)->add_block(residual, i0o, i0, i2);
            }
        break;
      case Excitations::airs:
        for (auto& i1 : virt_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i0, i1)) continue;
              unique_ptr<double[]> residual = rtensor->get_block(i0o, i0, i1);
              unique_ptr<double[]> denom    = dtensor->get_block(i0o, i0, i1);
              const size_t blocksize = rtensor->get_size(i0o, i0, i1);
              if (imag_) {
                for (size_t j = 0; j != blocksize; ++j) {
                  residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
                }
              } else {
                for (size_t j = 0; j != blocksize; ++j) {
                  residual[j] *= -(1.0 / (denom[j] + shift_));
                }
              }
              data_[istate]->at(iext)->add_block(residual, i0o, i0, i1);
            }
        break;
      case Excitations::arst:
        for (auto& i1 : virt_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i1)) continue;
            unique_ptr<double[]> residual = rtensor->get_block(i0o, i1);
            unique_ptr<double[]> denom    = dtensor->get_block(i0o, i1);
            const size_t blocksize = rtensor->get_size(i0o, i1);
            if (imag_) {
              for (size_t j = 0; j != blocksize; ++j) {
                residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
              }
            } else {
              for (size_t j = 0; j != blocksize; ++j) {
                residual[j] *= -(1.0 / (denom[j] + shift_));
              }
            }
            data_[istate]->at(iext)->add_block(residual, i0o, i1);
          }
        break;
      case Excitations::rist:
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i0)) continue;
            unique_ptr<double[]> residual = rtensor->get_block(i0o, i0);
            unique_ptr<double[]> denom    = dtensor->get_block(i0o, i0);
            const size_t blocksize = rtensor->get_size(i0o, i0);
            if (imag_) {
              for (size_t j = 0; j != blocksize; ++j) {
                residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
              }
            } else {
              for (size_t j = 0; j != blocksize; ++j) {
                residual[j] *= -(1.0 / (denom[j] + shift_));
              }
            }
            data_[istate]->at(iext)->add_block(residual, i0o, i0);
          }
        break;
      case Excitations::aibj:
        for (int ist = 0; ist != nstates_; ++ist) {
          if (!sssr_ || ist == istate) {
            const int pos = iext + (sssr_ ? 0 : ist);
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!denom_[istate]->at(pos)->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<double[]> residual = r->data(istate)->at(pos)->get_block(i0, i1, i2, i3);
                    unique_ptr<double[]> denom    = denom_[istate]->at(pos)->get_block(i0, i1, i2, i3);
                    const size_t blocksize = r->data(istate)->at(pos)->get_size(i0, i1, i2, i3);
                    if (imag_) {
                      for (size_t j = 0; j != blocksize; ++j) {
                        residual[j] *= -(1.0 / (denom[j] * denom[j] + shift2));
                      }
                    } else {
                      for (size_t j = 0; j != blocksize; ++j) {
                        residual[j] *= -(1.0 / (denom[j] + shift_));
                      }
                    }
                    data_[istate]->at(pos)->add_block(residual, i0, i1, i2, i3);
                  }
          }
        }
        break;
    }
  }
  mpi__->barrier();
}


void Orthogonal_Basis::print_convergence(shared_ptr<const Orthogonal_Basis> s, shared_ptr<const Orthogonal_Basis> r, const int istate, const int iter, const double error, const double timing) const {
  shared_ptr<MultiTensor_<double>> tcovar = get_contravariant(istate, true);
  double etot = 0.0;
  vector<double> esector(Excitations::total);

  cout << setw(8) << iter;
  for (int ist = 0; ist != nstates_; ++ist) {
    if (!sssr_ || istate == ist) {
      const int pos = Excitations::aibj + (sssr_ ? 0 : ist);
      esector[Excitations::aibj] += tcovar->at(pos)->dot_product(s->data(istate)->at(pos)) + tcovar->at(pos)->dot_product(r->data(istate)->at(pos));
    }
  }
  etot += esector[Excitations::aibj];
  cout << setprecision(6) << setw(12) << esector[Excitations::aibj];
  for (int iext = Excitations::arbs; iext != Excitations::aibj; ++iext) {
    esector[iext] = tcovar->at(iext)->dot_product(s->data(istate)->at(iext)) + tcovar->at(iext)->dot_product(r->data(istate)->at(iext));
    cout << setprecision(6) << setw(12) << esector[iext];
    etot += esector[iext];
  }

  cout << setw(15) << setprecision(8) << etot << setw(15) << error << setw(7) << setprecision(2) << timing << endl;
}


shared_ptr<VectorB> Orthogonal_Basis::vectorb(const int istate, const int iext) const {
  const size_t outsize = size(iext) * ((iext == Excitations::aibj && !sssr_) ? nstates_ : 1);
  auto out = make_shared<VectorB>(outsize);

  const size_t interm_size = (iext == Excitations::aibj) ? 0 : shalf_[iext]->ndim();
  const shared_ptr<Tensor_<double>> dtensor = data_[istate]->at(iext);
  switch(iext) {
    case Excitations::arbs:
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i1, i3)) continue;
            unique_ptr<double[]> data0 = dtensor->get_block(i0o, i1, i3);
            size_t iall = 0;
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
                  const size_t jall = j0o + interm_size * ((j1 - nocc_) + nvirt_ * (j3 - nocc_));
                  (*out)[jall] = data0[iall];
                }
          }
      break;
    case Excitations::arbi:
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
              unique_ptr<double[]> data0 = dtensor->get_block(i0o, i1, i2, i3);
              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
                      const size_t jall = j0o + interm_size * ((j1 - nocc_) + nvirt_ * ((j2 - ncore_) + nclo_ * (j3 - nocc_)));
                      (*out)[jall] = data0[iall];
                    }
            }
      break;
    case Excitations::airj:
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[iext]) {
              if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
              unique_ptr<double[]> data0 = dtensor->get_block(i0o, i0, i1, i2);
              size_t iall = 0;
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                    for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
                      const size_t jall = (j0 - ncore_) + nclo_ * ((j1 - nocc_) + nvirt_ * ((j2 - ncore_) + nclo_ * j0o));
                      (*out)[jall] = data0[iall];
                    }
            }
      break;
    case Excitations::risj:
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i0, i2)) continue;
            unique_ptr<double[]> data0 = dtensor->get_block(i0o, i0, i2);
            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
                  const size_t jall = j0 - ncore_ + nclo_ * ((j2 - ncore_) + nclo_ * j0o);
                  (*out)[jall] = data0[iall];
                }
          }
      break;
    case Excitations::airs:
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[iext]) {
            if (!dtensor->is_local(i0o, i0, i1)) continue;
            unique_ptr<double[]> data0 = dtensor->get_block(i0o, i0, i1);
            size_t iall = 0;
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
                for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
                  const size_t jall = j0 - ncore_ + nclo_ * (j1 - nocc_ + nvirt_ * j0o);
                  (*out)[jall] = data0[iall];
                }
          }
      break;
    case Excitations::arst:
      for (auto& i1 : virt_)
        for (auto& i0o : interm_[iext]) {
          if (!dtensor->is_local(i0o, i1)) continue;
          unique_ptr<double[]> data0 = dtensor->get_block(i0o, i1);
          size_t iall = 0;
          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
            for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
              const size_t jall = j1 - nocc_ + nvirt_ * j0o;
              (*out)[jall] = data0[iall];
            }
        }
      break;
    case Excitations::rist:
      for (auto& i0 : closed_)
        for (auto& i0o : interm_[iext]) {
          if (!dtensor->is_local(i0o, i0)) continue;
          unique_ptr<double[]> data0 = dtensor->get_block(i0o, i0);
          size_t iall = 0;
          for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0)
            for (int j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++iall) {
              const size_t jall = j0 - ncore_ + nclo_ * j0o;
              (*out)[jall] = data0[iall];
            }
        }
      break;
    case Excitations::aibj:
      int ioffset = 0;
      for (int ist = 0; ist != nstates_; ++ist) {
        if (!sssr_ || ist == istate) {
          const int pos = iext + (sssr_ ? 0 : ist);
           for (auto& i3 : virt_)
             for (auto& i2 : closed_)
               for (auto& i1 : virt_)
                 for (auto& i0 : closed_) {
                   if (!dtensor->is_local(i0, i1, i2, i3)) continue;
                   unique_ptr<double[]> data0 = data_[istate]->at(pos)->get_block(i0, i1, i2, i3);
                   size_t iall = 0;
                   for (int j3 = i3.offset()-nocc_; j3 != i3.offset()+i3.size()-nocc_; ++j3)
                     for (int j2 = i2.offset()-ncore_; j2 != i2.offset()+i2.size()-ncore_; ++j2)
                       for (int j1 = i1.offset()-nocc_; j1 != i1.offset()+i1.size()-nocc_; ++j1)
                         for (int j0 = i0.offset()-ncore_; j0 != i0.offset()+i0.size()-ncore_; ++j0, ++iall) {
                           const size_t jall = j0 + nclo_ * (j1 + nvirt_ * (j2 + nclo_ * j3));
                           (*out)[jall + ioffset] = data0[iall];
                         }
                 }
           ioffset += size(iext);
        }
      }
      break;
  }
  out->allreduce();
  return out;
}


// copy of the functions in SpinFreeMethod to initialize transformed tensors
void Orthogonal_Basis::loop_over(function<void(const Index&, const Index&, const Index&, const Index&)> func) const {
  for (auto& i3 : virt_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i2 : active_)
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i0 : active_)
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_) {
          func(i0, i1, i2, i3);
          func(i2, i1, i0, i3);
        }
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i0 : active_)
        for (auto& i2 : closed_)
          func(i0, i1, i2, i3);
}


shared_ptr<Tensor_<double>> Orthogonal_Basis::init_amplitude() const {
  unordered_set<size_t> sparse;
  auto put = [&sparse](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    sparse.insert(generate_hash_key(i0, i1, i2, i3));
  };
  loop_over(put);

  IndexRange occ(closed_);   occ.merge(active_);
  IndexRange virt(active_);  virt.merge(virt_);

  return make_shared<Tensor_<double>>(vector<IndexRange>{occ, virt, occ, virt}, /*kramers*/false, sparse, /*alloc*/true);
}


#endif
