//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: scf_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/scf_base.h>
#include <src/wfn/zreference.h>
#include <src/util/timer.h>
#include <src/util/math/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace bagel;


template <typename MatType, typename OvlType, typename HcType, class Enable>
SCF_base_<MatType, OvlType, HcType, Enable>::SCF_base_(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> re, const bool need_schwarz)
 : Method(idat, geom, re), eig_(geom->nbasis()) {

  // if this is called by Opt
  do_grad_ = idata_->get<bool>("_gradient", false);
  // enable restart capability
  restart_ = idata_->get<bool>("restart", false);
  dofmm_   = !(geom_->fmm() == nullptr);

  Timer scfb;
  overlap_ = make_shared<const OvlType>(geom_);
  scfb.tick_print("Overlap matrix");

  hcore_ = make_shared<HcType>(geom_, geom_->hcoreinfo());
  scfb.tick_print("Hcore matrix");

  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata_->get<int>("diis_start", 1);
  diis_size_ = idata_->get<int>("diis_size", 5);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh_scf", thresh_scf_);

  if (dofmm_) {
    fmm_ = make_shared<const FMM>(idata_, geom);
    const bool fmmk = idata_->get<bool>("FMM-K", false);
    if (fmmk)
      fmmK_ = make_shared<const FMM>(idata_, geom, true);
  }
  multipole_print_ = idata_->get<int>("multipole", 1);
  dma_print_ = idata_->get<int>("dma", 0);

  const int ncharge = idata_->get<int>("charge", 0);
  const int nopen   = idata_->get<int>("nopen", (geom_->nele()-ncharge)%2);
  nocc_ = (geom_->nele()-ncharge+nopen)/2;
  noccB_ = nocc_ - nopen;

  if (nocc_+noccB_ != geom_->nele()-ncharge) throw runtime_error("nocc and nopen are not consistently specified");

  tildex_ = overlap_->tildex(thresh_overlap_);

  scfb.tick_print("Overlap orthog");

  if (need_schwarz) {
    init_schwarz();
    scfb.tick_print("Schwarz matrix");
  }

  // if ref is passed to this
  if (re != nullptr) {
    get_coeff(re);
  }
  cout << endl;
}


template <typename MatType, typename OvlType, typename HcType, class Enable>
void SCF_base_<MatType, OvlType, HcType, Enable>::init_schwarz() {
  schwarz_ = geom_->schwarz();
}


// Specialized for GIAO
template <>
void SCF_base_<ZMatrix, ZOverlap, ZHcore, enable_if<true>::type>::get_coeff(const shared_ptr<const Reference> ref) {
  auto cref = dynamic_pointer_cast<const ZReference>(ref);
  assert(cref);
  coeff_ = cref->zcoeff();
}

template class bagel::SCF_base_<Matrix, Overlap, Hcore>;
template class bagel::SCF_base_<ZMatrix, ZOverlap, ZHcore>;

BOOST_CLASS_EXPORT_IMPLEMENT(SCF_base)
BOOST_CLASS_EXPORT_IMPLEMENT(SCF_base_London)

