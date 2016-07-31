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

}


void DistFCI::compute_rdm12(const int ist, const int jst) {

}


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> DistFCI::rdm34(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>(); 
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> DistFCI::rdm12_alpha(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>();
}


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> DistFCI::rdm34_alpha(const int ist, const int jst) const {
  return tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>();
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
DistFCI::compute_rdm12_from_civec(shared_ptr<const DistCivec>, shared_ptr<const DistCivec>) const {
  return tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>();
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
DistFCI::compute_rdm12_av_from_dvec(shared_ptr<const DistDvec>, shared_ptr<const DistDvec>, shared_ptr<const Determinants> o) const {
  return tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>();
}


shared_ptr<Dvec> DistFCI::rdm1deriv(const int istate) const {
  return nullptr;
}


shared_ptr<Dvec> DistFCI::rdm2deriv(const int istate) const {
  return nullptr;
}


tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>
DistFCI::rdm34deriv(const int istate, shared_ptr<const Matrix> fock, const size_t offset, const size_t size) const {
  return tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>();
}


pair<shared_ptr<Matrix>, VectorB> DistFCI::natorb_convert() {
  return pair<shared_ptr<Matrix>, VectorB>();
}
