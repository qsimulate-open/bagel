//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/multisite/multisite.h>
#include <src/mat1e/overlap.h>
#include <src/util/muffle.h>

using namespace std;
using namespace bagel;

MultiSite::MultiSite(shared_ptr<const PTree> input, shared_ptr<const Reference> ref, const int nsites) : input_(input), nsites_(nsites) {

}

void MultiSite::compute() {

}


