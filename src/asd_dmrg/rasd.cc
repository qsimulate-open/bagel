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
#include <src/ras/form_sigma.h>

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
    //cout << "RAS[" << nras[0] << "," << nras[1] << "," << nras[2] << "](" << input->get<int>("max_holes") << "h" << input->get<int>("max_particles") << "p)" << endl;
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

shared_ptr<Matrix> RASD::compute_sigma2e(shared_ptr<const RASDvec> cc, shared_ptr<const MOFile> jop) const {
  const int nstates = cc->ij();
  // Maybe batchsize should be an attribute of RASD
  FormSigmaRAS form_2e(input_->get_child("ras")->get<int>("batchsize", 512));
  shared_ptr<const RASDvec> sigma = form_2e(cc, nullptr, jop->mo2e(), vector<int>(nstates, static_cast<int>(false)));

  auto out = make_shared<Matrix>(nstates, nstates);
  for (int i = 0; i < nstates; ++i) {
    for (int j = 0; j < i; ++j)
      out->element(i,j) = out->element(j,i) = cc->data(i)->dot_product(*sigma->data(j));
    out->element(i,i) = cc->data(i)->dot_product(*sigma->data(i));
  }
  return out;
}

shared_ptr<Matrix> RASD::compute_spin(shared_ptr<const RASDvec> cc) const {
  const int nstates = cc->ij();
  shared_ptr<const RASDvec> sigma = cc->spin();
  auto out = make_shared<Matrix>(nstates, nstates);
  for (int i = 0; i < nstates; ++i) {
    for (int j = 0; j < i; ++j)
      out->element(i,j) = out->element(j,i) = cc->data(i)->dot_product(*sigma->data(j));
    out->element(i,i) = cc->data(i)->dot_product(*sigma->data(i));
  }
  return out;
}

shared_ptr<DMRG_Block> RASD::compute_first_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref) {
  map<BlockKey, shared_ptr<const RASDvec>> states;
  map<BlockKey, shared_ptr<const Matrix>> h_2e;
  map<BlockKey, shared_ptr<const Matrix>> spinmap;
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
      shared_ptr<const Matrix> hamiltonian_2e = compute_sigma2e(civecs, ras->jop());
      shared_ptr<const Matrix> spinmatrix = compute_spin(civecs);

      // Combines data for vectors with the same nelea and neleb
      auto organize_data = [&states, &h_2e, &spinmap] (shared_ptr<const RASDvec> civecs, shared_ptr<const Matrix> ham2e, shared_ptr<const Matrix> spinmat) {
        BlockKey key(civecs->det()->nelea(), civecs->det()->neleb());
        if (states.find(key) == states.end()) {
          assert(h_2e.find(key) == h_2e.end());
          states.emplace(key, civecs);
          h_2e.emplace(key, ham2e);
          spinmap.emplace(key, spinmat);
        }
        else {
          assert(h_2e.find(key) != h_2e.end());
          vector<shared_ptr<RASCivec>> tmpvecs = states[key]->dvec();
          vector<shared_ptr<RASCivec>> new_vecs = civecs->dvec();
          tmpvecs.insert(tmpvecs.end(), new_vecs.begin(), new_vecs.end());
          states[key] = make_shared<RASDvec>(tmpvecs);

          shared_ptr<Matrix> tmp2e = h_2e[key]->copy();
          const int oldsize = tmp2e->ndim();
          const int newsize = ham2e->ndim();
          tmp2e->resize(oldsize + newsize, oldsize + newsize);
          tmp2e->copy_block(oldsize, oldsize, newsize, newsize, *ham2e);
          h_2e[key] = tmp2e;

          shared_ptr<Matrix> tmpspin = spinmap[key]->copy();
          tmpspin->resize(oldsize + newsize, oldsize + newsize);
          tmpspin->copy_block(oldsize, oldsize, newsize, newsize, *spinmat);
          spinmap[key] = tmpspin;
        }
      };

      organize_data(civecs, hamiltonian_2e, spinmatrix);

      for (int i = 0; i < spin; ++i) {
        shared_ptr<RASDvec> tmpvec = civecs->spin_lower();
        for (auto& vec : tmpvec->dvec())
          vec->normalize();

        organize_data(tmpvec, hamiltonian_2e, spinmatrix);
        civecs = tmpvec;
      }
    }
    const int nstates = inp->get<int>("nstate");
    cout << "      - charge: " << charge << ", nspin: " << spin << ", nstates: " << nstates
                                    << fixed << setw(10) << setprecision(2) << rastime.tick() << endl;
  }

  GammaForestASD<RASDvec> forest(states);
  return make_shared<DMRG_Block>(move(forest), h_2e, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact()));
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
      map<BlockKey, shared_ptr<const RASDvec>> civecs = diagonalize_site_RDM(prod_ras->civectors());
      map<BlockKey, shared_ptr<const Matrix>> h_2e;
      map<BlockKey, shared_ptr<const Matrix>> spinmap;
      //for_each(civecs.begin(), civecs.end(), [] (pair<BlockKey, shared_ptr<const RASDvec>> c) { c.second->print(); });
      for (auto& ici : civecs) {
        h_2e.emplace(ici.first, compute_sigma2e(ici.second, prod_ras->jop()->monomer_jop<0>()));
        spinmap.emplace(ici.first, compute_spin(ici.second));
      }

      for (int i = 0; i < nstates_; ++i)
        sweep_energies_[i].push_back(prod_ras->energy(i));

      GammaForestASD<RASDvec> forest(civecs);
      return make_shared<DMRG_Block>(move(forest), h_2e, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact()));
    }
    else {
      throw logic_error("Full DMRG sweep not yet implemented!");
      return nullptr;
    }
  }
}

map<BlockKey, shared_ptr<const RASDvec>> RASD::diagonalize_site_RDM(const vector<shared_ptr<ProductRASCivec>>& civecs) const {
  assert(civecs.size()==nstates_);

  // store the coefficients and sector bases by block
  map<BlockKey, pair<shared_ptr<Matrix>, vector<pair<int, shared_ptr<const RASBlockVectors>>>>> basisdata;

  // arrange all the 'singular values' to get the best ones
  multimap<double, pair<BlockKey, const int>> singular_values;

  // non-orthogonal basis to do the diagonalization
  map<BlockKey, vector<shared_ptr<const RASBlockVectors>>> cibasis;
  for (auto& ivec : civecs)
    for (auto& isec : ivec->sectors())
      cibasis[isec.first].push_back(isec.second);

  // RDM is block diagonal by sector
  for (auto& isec : cibasis) {
    // form basis for sector
    vector<pair<int, shared_ptr<const RASBlockVectors>>> sectorbasis;
    {
      int offset = 0;
      for (auto& i : isec.second) {
        sectorbasis.emplace_back(offset, i);
        offset += i->mdim();
      }
    }
    const size_t bsize = accumulate(isec.second.begin(), isec.second.end(), 0ul, [] (size_t a, shared_ptr<const Matrix> m) { return a+m->mdim(); });

    // build overlap matrix
    Matrix overlap(bsize, bsize);
    for (auto i = sectorbasis.begin(); i != sectorbasis.end(); ++i) {
      shared_ptr<const Matrix> imat = i->second;
      for (auto j = sectorbasis.begin(); j != i; ++j) {
        if (i->second->left_state().key()==j->second->left_state().key()) {
          shared_ptr<const Matrix> jmat = j->second;
          assert(imat->ndim()==jmat->ndim());
          dgemm_("T", "N", imat->mdim(), jmat->mdim(), imat->ndim(), 1.0, imat->data(), imat->ndim(), jmat->data(), jmat->ndim(),
                                                                     0.0, overlap.element_ptr(i->first, j->first), overlap.ndim());
        }
      }
      dgemm_("T", "N", imat->mdim(), imat->mdim(), imat->ndim(), 1.0, imat->data(), imat->ndim(), imat->data(), imat->ndim(),
                                                                 0.0, overlap.element_ptr(i->first, i->first), overlap.ndim());
    }
    overlap.fill_upper();

    // build rdm
    Matrix rdm(bsize, bsize);
    for (int ist = 0; ist < nstates_; ++ist) {
      shared_ptr<const Matrix> sector = civecs[ist]->sector(isec.first);
      Matrix tmp(sector->mdim(), bsize);
      for (auto& ib : sectorbasis)
        dgemm_("T", "N", sector->mdim(), ib.second->mdim(), sector->ndim(), 1.0, sector->data(), sector->ndim(), ib.second->data(), ib.second->ndim(),
                                                                            0.0, tmp.element_ptr(0, ib.first), tmp.ndim());

      dgemm_("T", "N", bsize, bsize, sector->mdim(), weights_[ist], tmp.data(), tmp.ndim(), tmp.data(), tmp.ndim(),
                                                               1.0, rdm.data(), rdm.ndim());
    }

    Matrix orthonormalize(*overlap.tildex(1.0e-12));
    auto best_states = make_shared<Matrix>(orthonormalize % rdm * orthonormalize);
    VectorB eigs(best_states->ndim());
    best_states->diagonalize(eigs);
    best_states = make_shared<Matrix>(orthonormalize * *best_states);
    basisdata.emplace(isec.first, make_pair(best_states, sectorbasis));
    for (int i = 0; i < eigs.size(); ++i)
      singular_values.emplace(eigs(i), make_pair(isec.first, i));
  }

  // pre-organize output
  map<BlockKey, vector<shared_ptr<RASCivec>>> output_vectors;

  // Pick the top singular values
  int nvectors = 0;
  double partial_trace = 0.0;
  for (auto i = singular_values.rbegin(); i != singular_values.rend(); ++i, ++nvectors) {
    if (nvectors==ntrunc_) break;
    BlockKey bk = i->second.first;
    const int position = i->second.second;

    shared_ptr<const Matrix> coeff = basisdata[bk].first;
    vector<pair<int, shared_ptr<const RASBlockVectors>>>& sectorbasis = basisdata[bk].second;

    shared_ptr<const RASDeterminants> det = sectorbasis.front().second->det();

    auto tmp = make_shared<RASCivec>(det);
    for (auto& ic : sectorbasis) {
      shared_ptr<const RASBlockVectors> vecs = ic.second;
      dgemv_("N", det->size(), vecs->mdim(), 1.0, vecs->data(), vecs->ndim(), coeff->element_ptr(ic.first, position), 1,
                                             1.0, tmp->data(), 1);
    }

    cout << "spin expectation of vector: " << tmp->spin_expectation() << endl;
    output_vectors[i->second.first].push_back(tmp);
    partial_trace += i->first;
  }
  cout << "  discarded weights: " << setw(12) << setprecision(8) << scientific <<  1.0 - partial_trace << fixed << endl;

  // Process into Dvecs
  map<BlockKey, shared_ptr<const RASDvec>> out;
  for (auto& dvec : output_vectors)
    out.emplace(dvec.first, make_shared<RASDvec>(dvec.second));

  return out;
}
