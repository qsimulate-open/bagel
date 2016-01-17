//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_distras.cc
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

#include <src/asd/asd_distras.h>
#include <src/ci/ras/form_sigma.h>

using namespace std;
using namespace bagel;


shared_ptr<DistRASDvec> ASD_DistRAS::form_sigma(shared_ptr<const DistRASDvec> ccvec, shared_ptr<const MOFile> jop) const {
  vector<shared_ptr<RASCivec>> tmpvec;
  for (auto& i : ccvec->dvec()) tmpvec.push_back(make_shared<RASCivec>(i));
  auto dvec = make_shared<RASDvec>(tmpvec);

  FormSigmaRAS form;
  vector<int> conv(ccvec->ij(), static_cast<int>(false));
  shared_ptr<const RASDvec> sigmavec = form(dvec, jop, conv);

  return make_shared<DistRASDvec>(sigmavec);
}

// TODO function not yet written
shared_ptr<DistRASDvec> ASD_DistRAS::form_sigma_1e(shared_ptr<const DistRASDvec> ccvec, const double* modata) const {
  throw logic_error("ASD_DistRAS::form_sigma_1e function not yet written");
  return nullptr;
}

tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> ASD_DistRAS::compute_rdm12_monomer(shared_ptr<const DistRASDvec> civec, const int i) const {
  throw logic_error("not yet implemented in ASD_DistRAS");
  return make_tuple(nullptr,nullptr);
}
shared_ptr<DistRASDvec> ASD_DistRAS::contract_I(shared_ptr<const DistRASDvec> A, shared_ptr<Matrix> adiabats, int ioff, int nstA, int nstB, int kst) const {
  return nullptr;
}
shared_ptr<DistRASDvec> ASD_DistRAS::contract_J(shared_ptr<const DistRASDvec> A, shared_ptr<Matrix> adiabats, int ioff, int nstA, int nstB, int kst) const {
  return nullptr;
}
