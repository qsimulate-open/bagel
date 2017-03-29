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


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> DistFCI::rdm34(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>();
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


shared_ptr<Dvec> DistFCI::rdm1deriv(const int istate) const {
  return nullptr;
}


shared_ptr<Dvec> DistFCI::rdm2deriv(const int istate) const {
  return nullptr;
}


shared_ptr<Matrix> DistFCI::rdm2deriv_offset(const int istate, const size_t dsize, const size_t offset, const bool parallel) const {
  return nullptr;
}


tuple<shared_ptr<Matrix>,shared_ptr<Matrix>,shared_ptr<Matrix>>
DistFCI::rdm3deriv(const int istate, shared_ptr<const Matrix> fock, const size_t offset, const size_t size, shared_ptr<const Matrix> fock_ebra_in) const {
  return tuple<shared_ptr<Matrix>,shared_ptr<Matrix>,shared_ptr<Matrix>>();
}

tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>
DistFCI::rdm34deriv(const int istate, shared_ptr<const Matrix> fock, const size_t offset, const size_t size) const {
  return tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>();
}
