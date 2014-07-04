//
// BAGEL - Parallel electron correlation program.
// Filename: rasd.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd_dmrg/rasd.h>

using namespace std;
using namespace bagel;

RASD::RASD(const shared_ptr<const PTree> input, shared_ptr<Dimer> dimer) : ASD_DMRG(input, dimer) { }

shared_ptr<DMRG_Block> RASD::compute_first_block(shared_ptr<PTree> input, shared_ptr<const Reference> ref) {
  return nullptr;
}

shared_ptr<DMRG_Block> RASD::grow_block(shared_ptr<PTree> input, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block> left) {
  return nullptr;
}

shared_ptr<DMRG_Block> RASD::decimate_block(shared_ptr<PTree> input, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block> system, shared_ptr<DMRG_Block> environment) {
  for (int i = 0; i < nstates_; ++i) {
    sweep_energies_[i].push_back(0.0);
  }
  return nullptr;
}
