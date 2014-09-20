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
#include <src/asd_dmrg/form_sigma.h>
#include <src/ras/apply_operator.h>
#include <src/dimer/dimer_jop.h>
#include <src/util/muffle.h>
#include <src/ras/rasci.h>

#define DEBUG

using namespace std;
using namespace bagel;

RASD::RASD(const shared_ptr<const PTree> input, shared_ptr<MultiSite> multisite) : ASD_DMRG(input, multisite) { }

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

// TODO how am I going to ensure that I get the right formatting in Jop?
shared_ptr<Matrix> RASD::compute_sigma2e(const vector<shared_ptr<ProductRASCivec>>& cc, shared_ptr<const DimerJop> jop) const {
  const int nstates = cc.size();
  FormSigmaProdRAS form_2e(input_->get_child("ras")->get<int>("batchsize", 512));
  auto jop_2e = make_shared<DimerJop>(cc.front()->space()->norb(), cc.front()->left()->norb(), jop->mo1e()->clone(), jop->mo2e()->copy());
  shared_ptr<const BlockOperators> blockops = cc.front()->left()->compute_block_ops(jop_2e);
  vector<shared_ptr<ProductRASCivec>> sigma = form_2e(cc, blockops, jop_2e, vector<bool>(nstates, false));

  auto out = make_shared<Matrix>(nstates, nstates);
  for (int i = 0; i < nstates; ++i) {
    for (int j = 0; j < i; ++j)
      out->element(i,j) = out->element(j,i) = cc.at(i)->dot_product(*sigma.at(j));
    out->element(i,i) = cc.at(i)->dot_product(*sigma.at(i));
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

shared_ptr<Matrix> RASD::compute_spin(vector<shared_ptr<ProductRASCivec>> cc) const {
  const int nstates = cc.size();

  vector<shared_ptr<const ProductRASCivec>> spin_vec;
  for (auto& c : cc)
    spin_vec.push_back(c->spin());

  auto out = make_shared<Matrix>(nstates, nstates);
  for (int i = 0; i < nstates; ++i) {
    for (int j = 0; j < i; ++j)
      out->element(i,j) = out->element(j,i) = cc.at(i)->dot_product(*spin_vec.at(j));
    out->element(i,i) = cc.at(i)->dot_product(*spin_vec.at(i));
  }
  return out;
}

shared_ptr<DMRG_Block1> RASD::compute_first_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref) {
  map<BlockKey, shared_ptr<const RASDvec>> states;
  map<BlockKey, shared_ptr<const Matrix>> hmap;
  map<BlockKey, shared_ptr<const Matrix>> spinmap;
  Timer rastime;

  bool append = false;

  for (auto& inp : inputs) {
    // finish preparing the input
    inp->put("nclosed", ref->nclosed());
    read_restricted(inp, 0);
    const int spin = inp->get<int>("nspin");
    const int charge = inp->get<int>("charge");
    {
      Muffle hide_cout("asd_dmrg.log", append);
      append = true;
      // RAS calculations
      auto ras = make_shared<RASCI>(inp, ref->geom(), ref);
      ras->compute();
      shared_ptr<const RASDvec> civecs = ras->civectors();
      shared_ptr<const Matrix> hamiltonian_2e = compute_sigma2e(civecs, ras->jop());
      shared_ptr<const Matrix> spinmatrix = compute_spin(civecs);

      // Combines data for vectors with the same nelea and neleb
      auto organize_data = [&states, &hmap, &spinmap] (shared_ptr<const RASDvec> civecs, shared_ptr<const Matrix> ham2e, shared_ptr<const Matrix> spinmat) {
        BlockKey key(civecs->det()->nelea(), civecs->det()->neleb());
        if (states.find(key) == states.end()) {
          assert(hmap.find(key) == hmap.end());
          states.emplace(key, civecs);
          hmap.emplace(key, ham2e);
          spinmap.emplace(key, spinmat);
        }
        else {
          assert(hmap.find(key) != hmap.end());
          vector<shared_ptr<RASCivec>> tmpvecs = states[key]->dvec();
          vector<shared_ptr<RASCivec>> new_vecs = civecs->dvec();
          tmpvecs.insert(tmpvecs.end(), new_vecs.begin(), new_vecs.end());
          states[key] = make_shared<RASDvec>(tmpvecs);

          const int newsize = ham2e->ndim();
          const int oldsize = hmap[key]->ndim();
          shared_ptr<Matrix> tmp2e = hmap[key]->resize(oldsize+newsize, oldsize+newsize);
          tmp2e->copy_block(oldsize, oldsize, newsize, newsize, *ham2e);
          hmap[key] = tmp2e;

          shared_ptr<Matrix> tmpspin = spinmap[key]->resize(oldsize+newsize, oldsize+newsize);
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
  forest.compute();
  return make_shared<DMRG_Block1>(move(forest), hmap, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact()));
}

shared_ptr<DMRG_Block1> RASD::grow_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block1> left, const int site) {
  map<BlockKey, vector<shared_ptr<ProductRASCivec>>> states;
  map<BlockKey, shared_ptr<const Matrix>> hmap;
  map<BlockKey, shared_ptr<const Matrix>> spinmap;

  shared_ptr<const DimerJop> jop;

  Timer growtime;
  for (auto& inp : inputs) {
    // finish preparing the input
    const int charge = inp->get<int>("charge");
    const int spin = inp->get<int>("nspin");
    inp->put("nclosed", ref->nclosed());
    read_restricted(inp, site);
    {
      //Muffle hide_cout;
      // ProductRAS calculations
      auto prod_ras = make_shared<ProductRASCI>(inp, ref, left);
      prod_ras->compute();
      vector<shared_ptr<ProductRASCivec>> civecs = prod_ras->civectors();
      jop = prod_ras->jop();
      BlockKey key(prod_ras->nelea(), prod_ras->neleb());
      states[key].insert(states[key].end(), civecs.begin(), civecs.end());

      for (int i = 0; i < spin; ++i) {
        vector<shared_ptr<ProductRASCivec>> tmpvecs;
        for (auto& c : civecs) {
          shared_ptr<ProductRASCivec> out = c->spin_lower();
          out->normalize();
          tmpvecs.push_back(out);
        }

        key = BlockKey(key.nelea-1, key.neleb+1);
        states[key].insert(states[key].end(), tmpvecs.begin(), tmpvecs.end());
        civecs = tmpvecs;
      }
    }
    const int nstates = inp->get<int>("nstate");
    cout << "      - charge: " << charge << ", nspin: " << spin << ", nstates: " << nstates
                                    << fixed << setw(10) << setprecision(2) << growtime.tick() << endl;
  }


  assert(jop);
  for (auto& cc : states) {
    hmap.emplace(cc.first, compute_sigma2e(cc.second, jop));
    spinmap.emplace(cc.first, compute_spin(cc.second));
  }

  GammaForestProdASD forest(states);
  cout << "construct forest" << endl;
  forest.compute();
  cout << "compute GammaForestProdASD" << endl;

  shared_ptr<Matrix> coeff = ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact())->merge(states.begin()->second.front()->left()->coeff());
  return make_shared<DMRG_Block1>(move(forest), hmap, spinmap, coeff);
}

shared_ptr<DMRG_Block1> RASD::decimate_block(shared_ptr<PTree> input, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block1> system, shared_ptr<DMRG_Block1> environment, const int site) {
  Timer decimatetime;
  // assume the input is already fully formed, this may be revisited later
  input->put("nclosed", ref->nclosed());
  read_restricted(input, site);
  {
    Muffle hide_cout("asd_dmrg.log", true);
    // ProductRAS calculations
    if (!system) {
      auto prod_ras = make_shared<ProductRASCI>(input, ref, environment);
      prod_ras->compute();
      // diagonalize RDM to get RASCivecs
      map<BlockKey, shared_ptr<const RASDvec>> civecs = diagonalize_site_RDM(prod_ras->civectors(), perturb_);
      map<BlockKey, shared_ptr<const Matrix>> hmap;
      map<BlockKey, shared_ptr<const Matrix>> spinmap;
      for (auto& ici : civecs) {
        hmap.emplace(ici.first, compute_sigma2e(ici.second, prod_ras->jop()->monomer_jop<0>()));
        spinmap.emplace(ici.first, compute_spin(ici.second));
      }

      for (int i = 0; i < nstate_; ++i)
        sweep_energies_[i].push_back(prod_ras->energy(i));

      GammaForestASD<RASDvec> forest(civecs);
      forest.compute();
      return make_shared<DMRG_Block1>(move(forest), hmap, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact()));
    }
    else {
      auto block_pair = make_shared<DMRG_Block2>(system, environment);
#if 0
      auto prod_ras = make_shared<ProductRASCI>(input, ref, block_pair);
      prod_ras->compute();
      map<BlockKey, vector<shared_ptr<ProductRASCivec>>> civecs = diagonalize_site_and_block_RDM(prod_ras->civectors(), perturb_);

      const int norb = civecs.begin()->second.front()->det()->norb() + system->norb();
      const int nsysorb = civecs.begin()->second.front()->norb() + system->norb();
      shared_ptr<DimerJop> jop;
      { // make an appropriate DimerJop
        shared_ptr<const DimerJop> tmpjop = prod_ras->jop();
        auto mo1e = make_shared<CSymMatrix>(nsysorb);
        auto mo2e = make_shared<Matrix>(nsysorb*nsysorb, nsysorb*nsysorb);
        btas::TensorView mo2eView = btas::make_rview(btas::CRange<4>(norb,norb,norb,norb), tmpjop->mo2e()->storage());
        btas::TensorView mo2eSlice = btas::make_rview(mo2eView.range().slice({0, 0, 0, 0}, {nsysorb, nsysorb, nsysorb}), mo2eView.storage());
        copy(mo2eSlice.begin(), mo2eSlice.end(), mo2e->begin());
        jop = make_shared<DimerJop>(civecs.begin()->second.front()->det()->norb(), system->norb(), mo1e, mo2e);
      }

      map<BlockKey, shared_ptr<const Matrix>> hmap;
      map<BlockKey, shared_ptr<const Matrix>> spinmap;
      for (auto& c : civecs) {
        hmap.emplace(c.first, compute_sigma2e(civecs, jop));
        spinmap.emplace(c.first, compute_spin(civecs));
      }

      GammaForestProdASD forest(civecs);
      forest.compute();
      return make_shared<DMRG_Block1>(move(forest), hmap, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact())->merge(system->coeff()));
#else
      throw logic_error("Multi-site sweeps not yet implemented.");
      return nullptr;
#endif
    }
  }
}

map<BlockKey, shared_ptr<const RASDvec>> RASD::diagonalize_site_RDM(const vector<shared_ptr<ProductRASCivec>>& civecs, const double perturbation) const {
  assert(civecs.size()==nstate_);

  // store the coefficients and sector bases by block: BlockKey corresponds to the block's information
  map<BlockKey, shared_ptr<Matrix>> basisdata;

  // arrange all the 'singular values' to get the best ones
  multimap<double, tuple<BlockKey, const int>> singular_values;

  // arrange all the determinant objects for later use
  map<BlockKey, shared_ptr<const RASDeterminants>> detmap;

  // non-orthogonal basis to do the diagonalization
  map<BlockKey, shared_ptr<const Matrix>> cibasis;

  // list of all the vectors that get coupled together in the density matrix
  // rho = \sum_n \sum_{x,y in outer_products[n]} |x> w_n <y|
  map<BlockKey, vector<tuple<double, shared_ptr<const Matrix>>>> outer_products;

  // construct the non-orthogonal basis and the outer product vectors
  {
    // first, collect all of the vectors that belong in the state vectors
    for (int ist = 0; ist < nstate_; ++ist) {
      for (auto& isec : civecs[ist]->sectors()) {
        if (detmap.find(isec.first)==detmap.end())
          detmap[isec.first] = isec.second->det();
        outer_products[isec.first].emplace_back(weights_[ist], isec.second);
      }
    }

    // add in perturbative correction
    if (perturbation != 0.0) {
      for (int ist = 0; ist < nstate_; ++ist) {
        for (auto& isec : civecs[ist]->sectors()) {
          // "diagonal" perturbation: sum_{i,j} [ (i^dagger j)_alpha (i^dagger_j)_beta ]
          apply_perturbation(isec.second, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, detmap, weights_[ist]*perturbation, outer_products);
          apply_perturbation(isec.second, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}, detmap, weights_[ist]*perturbation, outer_products);

          // spinflip perturbations
          apply_perturbation(isec.second, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}, detmap, weights_[ist]*perturbation, outer_products);
          apply_perturbation(isec.second, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, detmap, weights_[ist]*perturbation, outer_products);

          // ET/HT perturbations
          apply_perturbation(isec.second, {GammaSQ::AnnihilateAlpha}, detmap, weights_[ist]*perturbation, outer_products);
          apply_perturbation(isec.second, {GammaSQ::AnnihilateBeta}, detmap, weights_[ist]*perturbation, outer_products);
          apply_perturbation(isec.second, {GammaSQ::CreateAlpha}, detmap, weights_[ist]*perturbation, outer_products);
          apply_perturbation(isec.second, {GammaSQ::CreateBeta}, detmap, weights_[ist]*perturbation, outer_products);
        }
      }

    }

    for (auto& basis_sector : outer_products) {
      const int bsize = accumulate(basis_sector.second.begin(), basis_sector.second.end(), 0, [] (int x, tuple<double, shared_ptr<const Matrix>> m) { return x+get<1>(m)->mdim(); });
      auto tmp_mat = make_shared<Matrix>(get<1>(basis_sector.second.front())->ndim(), bsize);
      int current = 0;
      for (auto& ib : basis_sector.second) {
        copy_n(get<1>(ib)->data(), get<1>(ib)->size(), tmp_mat->element_ptr(0, current));
        current += get<1>(ib)->mdim();
      }
      cibasis.emplace(basis_sector.first, tmp_mat);
    }
  }

  // RDM is block diagonal by sector
  for (auto& isec : cibasis) {
    shared_ptr<const Matrix> sectorbasis = isec.second;

    // build overlap matrix
    Matrix overlap(*sectorbasis % *sectorbasis);

    // build rdm
    Matrix rdm = *overlap.clone();
    for (auto& op : outer_products[isec.first]) {
      Matrix tmp(*get<1>(op) % *sectorbasis);
      dgemm_("T", "N", rdm.ndim(), rdm.mdim(), tmp.ndim(), get<0>(op), tmp.data(), tmp.ndim(), tmp.data(), tmp.ndim(), 1.0, rdm.data(), rdm.ndim());
    }

    Matrix orthonormalize(*overlap.tildex(1.0e-12));
    if (orthonormalize.mdim() > 0) {
      auto best_states = make_shared<Matrix>(orthonormalize % rdm * orthonormalize);
      VectorB eigs(best_states->ndim());
      best_states->diagonalize(eigs);
      best_states = make_shared<Matrix>(orthonormalize * *best_states);
      basisdata.emplace(isec.first, best_states);
      for (int i = 0; i < eigs.size(); ++i)
        singular_values.emplace(eigs(i), make_tuple(isec.first, i));
    }
  }

  // pre-organize output
  map<BlockKey, vector<shared_ptr<RASCivec>>> output_vectors;

  // Pick the top singular values
  int nvectors = 0;
  double partial_trace = 0.0;
  for (auto i = singular_values.rbegin(); i != singular_values.rend(); ++i, ++nvectors) {
    if (nvectors==ntrunc_) break;
    BlockKey bk = get<0>(i->second);
    const int position = get<1>(i->second);

    shared_ptr<const Matrix> coeff = basisdata[bk];
    shared_ptr<const Matrix> sectorbasis = cibasis[bk];
    shared_ptr<const RASDeterminants> det = detmap[bk];

    auto tmp = make_shared<RASCivec>(det);
    dgemv_("N", det->size(), sectorbasis->mdim(), 1.0, sectorbasis->data(), sectorbasis->ndim(), coeff->element_ptr(0, position), 1, 1.0, tmp->data(), 1);

    output_vectors[get<0>(i->second)].push_back(tmp);
    partial_trace += i->first;
  }
  cout << "  discarded weights: " << setw(12) << setprecision(8) << scientific <<  1.0 - partial_trace << fixed << endl;

  // Process into Dvecs: These vectors will become blocks so now BlockKey should be a descriptor of the block vector
  map<BlockKey, shared_ptr<const RASDvec>> out;
  for (auto& dvec : output_vectors) {
    auto tmp = make_shared<RASDvec>(dvec.second);
    out.emplace(BlockKey(tmp->det()->nelea(), tmp->det()->neleb()), tmp);
  }

  return out;
}

void RASD::apply_perturbation(shared_ptr<const RASBlockVectors> cc, vector<GammaSQ> oplist, map<BlockKey, shared_ptr<const RASDeterminants>>& detmap,
                        const double weight, map<BlockKey, vector<tuple<double, shared_ptr<const Matrix>>>>& outer_products) const
{
  pair<int, int> dele(0, 0);
  for (auto& op : oplist) {
    if (is_alpha(op))
      dele.first += (is_creation(op) ? 1 : -1);
    else
      dele.second += (is_creation(op) ? 1 : -1);
  }

  shared_ptr<const RASDeterminants> sdet = cc->det();
  const BlockKey basekey = cc->left_state().key();
  // block should have opposite action as dets
  const BlockKey Tkey(basekey.nelea - dele.first, basekey.neleb - dele.second);
  shared_ptr<const RASDeterminants> tdet;
  if (detmap.find(Tkey)!=detmap.end()) {
    tdet = detmap[Tkey];
  }
  else {
    if (sdet->nelea()+dele.first < 0 || sdet->neleb()+dele.second < 0) return;
    tdet = sdet->clone(sdet->nelea()+dele.first, sdet->neleb()+dele.second);
    detmap[Tkey] = tdet;
  }

  if (tdet->size()!=0) {
    ApplyOperator apply;
    const int nstates = cc->mdim();
    const int norb = tdet->norb();
    for (int ist = 0; ist < nstates; ++ist) {
      auto tmp = make_shared<Matrix>(tdet->size(), 1);
      RASCivecView view(tdet, tmp->data());
      if (oplist.size()==1) {
        for (int p = 0; p < norb; ++p)
          apply(1.0, cc->civec(ist), view, oplist, {p});
      }
      else if (oplist.size()==2) {
        for (int p = 0; p < norb; ++p)
          for (int q = 0; q < norb; ++q)
            apply(1.0, cc->civec(ist), view, oplist, {p, q});
      }
      else {
        assert(false);
      }
      outer_products[Tkey].emplace_back(weight, tmp);
    }
  }
}
