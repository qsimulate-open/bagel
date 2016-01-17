//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_distcas.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/asd_distcas.h>
#include <src/ci/fci/dist_form_sigma.h>

using namespace std;
using namespace bagel;

shared_ptr<DistDvec> ASD_DistCAS::form_sigma(shared_ptr<const DistDvec> ccvec, std::shared_ptr<const MOFile> jop) const {
  FormSigmaDistFCI form;

  return form(ccvec, jop);
}

shared_ptr<DistDvec> ASD_DistCAS::form_sigma_1e(shared_ptr<const DistDvec> ccvec, const double* modata) const {
  throw logic_error("1e properties not yet implemented in ASD_DistCAS");

  return nullptr;
}

tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> ASD_DistCAS::compute_rdm12_monomer(shared_ptr<const DistDvec> civec, const int i) const {
  throw logic_error("not yet implemented in ASD_DistCAS");
  return make_tuple(nullptr,nullptr);
}

shared_ptr<DistDvec> ASD_DistCAS::contract_I(shared_ptr<const DistDvec> A, shared_ptr<Matrix> adiabats, int ioff, int nstA, int nstB, int kst) const {
  return nullptr;
}
shared_ptr<DistDvec> ASD_DistCAS::contract_J(shared_ptr<const DistDvec> A, shared_ptr<Matrix> adiabats, int ioff, int nstA, int nstB, int kst) const {
  return nullptr;
}
