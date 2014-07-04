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

vector<shared_ptr<PTree>> ASD_DMRG::prepare_growing_input(const int site) const {
  vector<shared_ptr<PTree>> out;

  shared_ptr<PTree> base_input = input_->get_child_optional(input_->get<string>("method"));
  if (!base_input) base_input = make_shared<PTree>();

  base_input->erase("charge");
  base_input->erase("spin");
  base_input->erase("nstate");

  auto space = input_->get_child_optional("spaces");
  if (!space)
    throw runtime_error("spaces must be speciified");
  for (auto& s : *space) {
    auto tmp = make_shared<PTree>(*base_input);
    for (auto& siter : *s) {
      tmp->put("charge", siter->get<string>("charge"));
      tmp->put("spin", siter->get<string>("spin"));
      tmp->put("nstate", siter->get<string>("nstate"));
    }
    out.push_back(tmp);
  }

  if ( (out.size() != 1) && (out.size() != nsites_) )
    throw runtime_error("Must specify either \"space\" or one per site");

  return out;
}

shared_ptr<PTree> ASD_DMRG::prepare_sweeping_input(const int site) const {
  return nullptr;
}
