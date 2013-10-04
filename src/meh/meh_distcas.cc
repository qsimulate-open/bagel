//
// BAGEL - Parallel electron correlation program.
// Filename: meh_distcas.cc
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

#include <src/meh/meh_distcas.h>
#include <src/fci/dist_form_sigma.h>

using namespace std;
using namespace bagel;

shared_ptr<DistDvec> MEH_DistCAS::form_sigma(shared_ptr<const DistDvec> ccvec, std::shared_ptr<const MOFile> jop) const {
  FormSigmaDistFCI form;

  return form(ccvec, jop);
}

shared_ptr<DistDvec> MEH_DistCAS::form_sigma_1e(shared_ptr<const DistDvec> ccvec, const double* modata) const {
  throw logic_error("1e properties not yet implemented in MEH_DistCAS");

  return shared_ptr<DistDvec>();
}
