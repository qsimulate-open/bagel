//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: determinants_base.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/ci/ciutil/determinants_base.h>

template class bagel::Determinants_base<bagel::FCIString>;
template class bagel::Determinants_base<bagel::RASString>;

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Determinants_base<bagel::FCIString>)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Determinants_base<bagel::RASString>)
