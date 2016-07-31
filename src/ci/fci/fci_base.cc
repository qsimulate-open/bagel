//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_base.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/ci/fci/fci_base.h>

//BOOST_CLASS_EXPORT_IMPLEMENT(bagel::FCI_base<Civec,Dvec>)
//BOOST_CLASS_EXPORT_IMPLEMENT(bagel::FCI_base<DistCivec,DistDvec>)

using namespace bagel;

template class FCI_base<Civec,Dvec>;
template class FCI_base<DistCivec,DistDvec>;
