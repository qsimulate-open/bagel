//
// BAGEL - Parallel electron correlation program.
// Filename: asd_distras.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <src/asd/asd_distras.h>
#include <src/ras/form_sigma.h>

using namespace std;
using namespace bagel;

ASD_DistRAS::ASD_DistRAS(const shared_ptr<const PTree> input, shared_ptr<Dimer> dimer, shared_ptr<DimerDistRAS> cispace)
 : ASD<DistRASDvec>(input, dimer, cispace) {}


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
