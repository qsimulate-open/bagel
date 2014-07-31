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
#include <src/asd_dmrg/gamma_forest_asd.h>
#include <src/asd_dmrg/product_rasci.h>

#define DEBUG

using namespace std;
using namespace bagel;

RASD::RASD(const shared_ptr<const PTree> input, shared_ptr<Dimer> dimer) : ASD_DMRG(input, dimer) { }

void RASD::read_restricted(shared_ptr<PTree> input, const int site) const {
  auto restricted = input_->get_child("restricted");

  auto read = [&input] (const shared_ptr<const PTree> inp, int current) {
    array<int, 3> nras = inp->get_array<int, 3>("orbitals");
    input->put("max_holes", inp->get<string>("max_holes"));
    input->put("max_particles", inp->get<string>("max_particles"));

    input->erase("active");
    auto parent = std::make_shared<PTree>();
    for (int i = 0; i < 3; ++i) {
      auto tmp = std::make_shared<PTree>();
      const int norb = nras[i];
      for (int i = 0; i < norb; ++i, ++current)
        tmp->push_back(current+1);
      parent->push_back(tmp);
    }
    input->add_child("active", parent);
#ifdef DEBUG
    cout << "RAS[" << nras[0] << "," << nras[1] << "," << nras[2] << "](" << input->get<int>("max_holes") << "h" << input->get<int>("max_particles") << "p)" << endl;
#endif
  };

  if (restricted->size() == 1)
    read(*restricted->begin(), input->get<int>("nclosed"));
  else if (restricted->size() == nsites_) {
    auto iter = restricted->begin();
    advance(iter, site);
    read(*iter, input->get<int>("nclosed"));
  }
  else
    throw runtime_error("Must specify either one set of restrictions for all sites, or one set per site");
}

shared_ptr<DMRG_Block> RASD::compute_first_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref) {
  map<BlockKey, shared_ptr<const RASDvec>> states;
  map<BlockKey, shared_ptr<const Matrix>> h_2e;
  Timer rastime;

  for (auto& inp : inputs) {
    // finish preparing the input
    inp->put("nclosed", ref->nclosed());
    read_restricted(inp, 0);
    const int spin = inp->get<int>("nspin");
    const int charge = inp->get<int>("charge");
    {
      Muffle hide_cout;
      // RAS calculations
      auto ras = make_shared<RASCI>(inp, ref->geom(), ref);
      ras->compute();
      shared_ptr<const RASDvec> civecs = ras->civectors();
      shared_ptr<const Matrix> hamiltonian_2e = ras->compute_sigma2e();
      hamiltonian_2e->print();

      // Combines data for vectors with the same nelea and neleb
      auto organize_data = [&states, &h_2e] (shared_ptr<const RASDvec> civecs, shared_ptr<const Matrix> ham2e) {
        BlockKey key(civecs->det()->nelea(), civecs->det()->neleb());
        if (states.find(key) == states.end()) {
          assert(h_2e.find(key) == h_2e.end());
          states.emplace(key, civecs);
          h_2e.emplace(key, ham2e);
        }
        else {
          assert(h_2e.find(key) != h_2e.end());
          vector<shared_ptr<RASCivec>> tmpvecs = states[key]->dvec();
          vector<shared_ptr<RASCivec>> new_vecs = civecs->dvec();
          tmpvecs.insert(tmpvecs.end(), new_vecs.begin(), new_vecs.end());
          states[key] = std::make_shared<RASDvec>(tmpvecs);

          shared_ptr<Matrix> tmp2e = h_2e[key]->copy();
          const int oldsize = tmp2e->ndim();
          const int newsize = ham2e->ndim();
          tmp2e->resize(oldsize + newsize, oldsize + newsize);
          tmp2e->copy_block(oldsize, oldsize, newsize, newsize, *ham2e);
          h_2e[key] = tmp2e;
        }
      };

      organize_data(civecs, hamiltonian_2e);

      for (int i = 0; i < spin; ++i) {
        shared_ptr<RASDvec> tmpvec = civecs->spin_lower();
        for (auto& vec : tmpvec->dvec())
          vec->normalize();

        organize_data(tmpvec, hamiltonian_2e);
        civecs = tmpvec;
      }
    }
    const int nstates = inp->get<int>("nstate");
    cout << "      - charge: " << charge << ", nspin: " << spin << ", nstates: " << nstates
                                    << fixed << setw(10) << setprecision(2) << rastime.tick() << endl;
  }

  GammaForestASD<RASDvec> forest(states);
  return make_shared<DMRG_Block>(move(forest), h_2e, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact()));
}

shared_ptr<DMRG_Block> RASD::grow_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block> left, const int site) {
  map<BlockKey, vector<shared_ptr<const ProductRASCivec>>> states;
  map<BlockKey, shared_ptr<const Matrix>> h_2e;

  Timer growtime;
  for (auto& inp : inputs) {
    // finish preparing the input
    const int charge = inp->get<int>("charge");
    const int spin = inp->get<int>("spin");
    inp->put("nclosed", ref->nclosed());
    read_restricted(inp, site);
    {
      //Muffle hide_cout;
      // ProductRAS calculations
      auto prod_ras = make_shared<ProductRASCI>(inp, ref, left);
      prod_ras->compute();
      vector<shared_ptr<ProductRASCivec>> civecs = prod_ras->civectors();
      //shared_ptr<const Matrix> hamiltonian_2e = prod_ras->compute_sigma2e();
      auto hamiltonian_2e = make_shared<Matrix>(civecs.size(), civecs.size());
      //hamiltonian_2e->print();

      // Combines data for vectors with the same nelea and neleb
      auto organize_data = [&states, &h_2e, &prod_ras] (vector<shared_ptr<ProductRASCivec>> civecs, shared_ptr<const Matrix> ham2e) {
        BlockKey key(prod_ras->nelea(), prod_ras->neleb());
        states[key].insert(states[key].end(), civecs.begin(), civecs.end());
        if (h_2e.find(key) == h_2e.end()) {
          h_2e.emplace(key, ham2e);
        }
        else {
          shared_ptr<Matrix> tmp2e = h_2e[key]->copy();
          const int oldsize = tmp2e->ndim();
          const int newsize = ham2e->ndim();
          tmp2e->resize(oldsize + newsize, oldsize + newsize);
          tmp2e->copy_block(oldsize, oldsize, newsize, newsize, *ham2e);
          h_2e[key] = tmp2e;
        }
      };

      organize_data(civecs, hamiltonian_2e);

#if 0
      for (int i = 0; i < spin; ++i) {
        shared_ptr<RASDvec> tmpvec = civecs->spin_lower();
        for (auto& vec : tmpvec->dvec())
          vec->normalize();

        organize_data(tmpvec, hamiltonian_2e);
        civecs = tmpvec;
      }
#endif
    }
    const int nstates = inp->get<int>("nstate");
    cout << "      - charge: " << charge << ", nspin: " << spin << ", nstates: " << nstates
                                    << fixed << setw(10) << setprecision(2) << growtime.tick() << endl;
  }

#if 1
  return nullptr;
#else
  return make_shared<DMRG_Block>(/*stuff*/);
#endif
}

shared_ptr<DMRG_Block> RASD::decimate_block(shared_ptr<PTree> input, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block> system, shared_ptr<DMRG_Block> environment, const int site) {
  //map<BlockKey, vector<shared_ptr<const ProductRASCivec>>> states;
  //map<BlockKey, shared_ptr<const Matrix>> h_2e;

  Timer decimatetime;
  // assume the input is already fully formed, this may be revisited later
  input->put("nclosed", ref->nclosed());
  read_restricted(input, site);
  {
    //Muffle hide_cout;
    // ProductRAS calculations
    if (!system) {
      auto prod_ras = make_shared<ProductRASCI>(input, ref, environment);
      prod_ras->compute();
      // diagonalize RDM to get RASCivecs
      vector<shared_ptr<ProductRASCivec>> civecs = prod_ras->civectors();
      //shared_ptr<const Matrix> hamiltonian_2e = prod_ras->compute_sigma2e();
      //auto hamiltonian_2e = make_shared<Matrix>(civecs.size(), civecs.size());
      //hamiltonian_2e->print();
      for (int i = 0; i < nstates_; ++i)
        sweep_energies_[i].push_back(prod_ras->energy(i));
    }
    else {
      throw logic_error("Full DMRG sweep not yet implemented!");
    }
  }
#if 1
  return nullptr;
#else
  return make_shared<DMRG_Block>(/*stuff*/);
#endif
}
