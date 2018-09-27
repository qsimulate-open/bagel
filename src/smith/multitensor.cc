//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multitensor.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/smith/multitensor.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class bagel::SMITH::MultiTensor_<double>;
template class bagel::SMITH::MultiTensor_<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_CLASS_EXPORT_IMPLEMENT(MultiTensor_<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(MultiTensor_<complex<double>>)

#endif
