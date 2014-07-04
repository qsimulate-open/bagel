//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg_base.cc
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

#include <iostream>

#include <src/asd_dmrg/asd_dmrg.h>

using namespace std;
using namespace bagel;

ASD_DMRG::ASD_DMRG(shared_ptr<const PTree> input, shared_ptr<const Dimer> dimer) : input_(input), dimer_(dimer) {
  // Hardcoded to dimer for now
  nsites_ = 2;

  nstates_ = input_->get<int>("nstates", 1);
  thresh_ = input_->get<double>("thresh", 1.0e-8);
  maxiter_ = input_->get<int>("maxiter", 10);

  energies_.resize(nstates_);
  sweep_energies_.resize(nstates_);
}

string ASD_DMRG::print_progress(const int position, const string left_symbol, const string right_symbol) const {
  stringstream out;
  for (int i = 0; i < position; ++i) out << left_symbol << " ";
  out << "** ";
  for (int i = position+1; i < nsites_; ++i) out << right_symbol << " ";

  return out.str();
}

shared_ptr<PTree> ASD_DMRG::prepare_input(const int site, const bool growing) const {
  return nullptr;
}
