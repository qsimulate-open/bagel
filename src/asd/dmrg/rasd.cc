//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rasd.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/asd/dmrg/rasd.h>
#include <src/asd/dmrg/gamma_forest_asd.h>
#include <src/asd/dmrg/product_rasci.h>
#include <src/asd/dmrg/form_sigma.h>
#include <src/asd/dimer/dimer_jop.h>
#include <src/ci/ras/form_sigma.h>
#include <src/ci/ras/apply_operator.h>
#include <src/ci/ras/rasci.h>
#include <src/util/muffle.h>

using namespace std;
using namespace bagel;

RASD::RASD(const shared_ptr<const PTree> input, shared_ptr<const Reference> ref) : ASD_DMRG(input, ref) { }

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
#ifdef HAVE_MPI_H
  out->synchronize();
#endif
  return out;
}

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
#ifdef HAVE_MPI_H
  out->synchronize();
#endif
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
#ifdef HAVE_MPI_H
  out->synchronize();
#endif
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
#ifdef HAVE_MPI_H
  out->synchronize();
#endif
  return out;
}

shared_ptr<DMRG_Block1> RASD::compute_first_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref) {
  map<BlockKey, shared_ptr<const RASDvec>> states;
  map<BlockKey, shared_ptr<const Matrix>> hmap;
  map<BlockKey, shared_ptr<const Matrix>> spinmap;
  Timer rastime;

  bool append = false;

  for (auto& inp : inputs) {
    const int spin = inp->get<int>("nspin");
    const int charge = inp->get<int>("charge");
    { // prepare the input
      inp->put("nclosed", ref->nclosed());
      inp->put("extern_nactele", true);
      inp->put("nactele", active_electrons_.at(0));
      read_restricted(inp, 0);
    }
    { // RAS calculations
      Muffle hide_cout("asd_dmrg.log", append);
      append = true;
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
        for (auto& vec : tmpvec->dvec()) {
          vec->normalize();
#ifdef HAVE_MPI_H
          vec->synchronize();
#endif
        }

        organize_data(tmpvec, hamiltonian_2e, spinmatrix);
        civecs = tmpvec;
      }
    }
    const int nstates = inp->get<int>("nstate");
    cout << "      - charge: " << charge << ", nspin: " << spin << ", nstates: " << nstates
                                    << fixed << setw(10) << setprecision(2) << rastime.tick() << endl;
  }

  GammaForestASD<RASDvec> forest(states);
  rastime.tick_print("construct forest");
  forest.compute();
  rastime.tick_print("compute forest");
  auto coeff = ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed() + ref->nact());
  rastime.tick_print("coeff ");
  auto out = make_shared<DMRG_Block1>(move(forest), hmap, spinmap, coeff);
  rastime.tick_print("construct dmrg");
  return out;
}

shared_ptr<DMRG_Block1> RASD::grow_block(vector<shared_ptr<PTree>> inputs, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block1> left, const int site) {
  map<BlockKey, vector<shared_ptr<ProductRASCivec>>> states;
  map<BlockKey, shared_ptr<const Matrix>> hmap;
  map<BlockKey, shared_ptr<const Matrix>> spinmap;

  shared_ptr<const DimerJop> jop;

  Timer growtime(2);
  for (auto& inp : inputs) {
    const int charge = inp->get<int>("charge");
    const int spin = inp->get<int>("nspin");
    { // prepare input
      inp->put("nclosed", ref->nclosed());
      inp->put("extern_nactele", true);
      const int nactele = accumulate(active_electrons_.begin(), active_electrons_.begin()+site+1, 0);
      inp->put("nactele", nactele);
      read_restricted(inp, site);
    }
    { // ProductRAS calculations
      Muffle hide_cout("asd_dmrg.log", true);
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
  map<BlockKey, vector<shared_ptr<ProductRASCivec>>> ortho_states;
  for (auto& cc : states) {
    // orthogonalize the civecs first just in case
    Matrix overlap(cc.second.size(), cc.second.size());
    for (int i = 0; i < cc.second.size(); ++i) {
      for (int j = 0; j < i; ++j)
        overlap(i,j) = overlap(j,i) = cc.second[i]->dot_product(*cc.second[j]);
      overlap(i,i) = cc.second[i]->dot_product(*cc.second[i]);
    }
#ifdef HAVE_MPI_H
    overlap.synchronize();
#endif
    Matrix ortho = *overlap.tildex(1e-11);

    vector<shared_ptr<ProductRASCivec>> tmpvec;

    for (int i = 0; i < ortho.mdim(); ++i) {
      auto tmp = cc.second.front()->clone();
      for (int j = 0; j < ortho.ndim(); ++j)
        tmp->ax_plus_y(ortho(j,i), *cc.second[j]);
#ifdef HAVE_MPI_H
      tmp->synchronize();
#endif
      tmpvec.push_back(tmp);
    }

    hmap.emplace(cc.first, compute_sigma2e(tmpvec, jop));
    spinmap.emplace(cc.first, compute_spin(tmpvec));
    ortho_states.emplace(cc.first, move(tmpvec));
  }
  growtime.tick_print("orthonormalize and collect individual states");

  GammaForestProdASD forest(ortho_states);
  growtime.tick_print("construct GammaForestProdASD");
  forest.compute();
  growtime.tick_print("compute forest");

  shared_ptr<Matrix> coeff = ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact())->merge(left->coeff());
  auto out = make_shared<DMRG_Block1>(move(forest), hmap, spinmap, coeff);
  growtime.tick_print("dmrg block");

  return out;
}

shared_ptr<DMRG_Block1> RASD::decimate_block(shared_ptr<PTree> input, shared_ptr<const Reference> ref, shared_ptr<DMRG_Block1> system, shared_ptr<DMRG_Block1> environment, const int site) {
  Timer decimatetime(2);
  { // prepare input
    input->put("nclosed", ref->nclosed());
    input->put("extern_nactele", true);
    const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), input->get<int>("charge"));
    input->put("nactele", nactele);
    read_restricted(input, site);
  }
  { // ProductRAS calculations
    Muffle hide_cout("asd_dmrg.log", true);
    if (!system) {
      auto prod_ras = make_shared<ProductRASCI>(input, ref, environment);
      prod_ras->compute();
      decimatetime.tick_print("ProductRASCI calculation");

      // diagonalize RDM to get RASCivecs
      map<BlockKey, shared_ptr<const RASDvec>> civecs = diagonalize_site_RDM(prod_ras->civectors(), perturb_);
      decimatetime.tick_print("diagonalize site RDM");

      map<BlockKey, shared_ptr<const Matrix>> hmap;
      map<BlockKey, shared_ptr<const Matrix>> spinmap;
      for (auto& ici : civecs) {
        hmap.emplace(ici.first, compute_sigma2e(ici.second, prod_ras->jop()->monomer_jop<0>()));
        spinmap.emplace(ici.first, compute_spin(ici.second));
      }
      decimatetime.tick_print("compute renormalized 2e energy and spin");

      for (int i = 0; i < nstate_; ++i)
        sweep_energies_[i].push_back(prod_ras->energy(i));

      GammaForestASD<RASDvec> forest(civecs);
      decimatetime.tick_print("construct GammaForestASD");

      forest.compute();
      decimatetime.tick_print("compute forest");

      auto out = make_shared<DMRG_Block1>(move(forest), hmap, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact()));
      decimatetime.tick_print("dmrg block");

      return out;
    }
    else {
      auto block_pair = make_shared<DMRG_Block2>(system, environment);
      decimatetime.tick_print("Build double block");
      auto prod_ras = make_shared<ProductRASCI>(input, ref, block_pair);
      prod_ras->compute();
      decimatetime.tick_print("ProductRASCI calculation");

      for (int i = 0; i < nstate_; ++i)
        sweep_energies_[i].push_back(prod_ras->energy(i));
      decimatetime.tick_print("add results to vector");

      map<BlockKey, vector<shared_ptr<ProductRASCivec>>> civecs = diagonalize_site_and_block_RDM(prod_ras->civectors(), perturb_);
      decimatetime.tick_print("diagonalize system RDM");

      const int nrasorb = civecs.begin()->second.front()->space()->norb();
      const int nsysorb = nrasorb + system->norb();
      const int norb = nsysorb + environment->norb();

      // make an appropriate DimerJop
      auto mo2e = make_shared<Matrix>(nsysorb*nsysorb, nsysorb*nsysorb);
      const btas::TensorView4<double> mo2eView = btas::make_view(btas::CRange<4>(norb,norb,norb,norb), prod_ras->jop()->mo2e()->storage());
      auto low = {0, 0, 0, 0};
      auto up = {nsysorb, nsysorb, nsysorb, nsysorb};
      const btas::TensorView4<double> mo2eSlice = btas::make_view(mo2eView.range().slice(low, up), mo2eView.storage());
      copy(mo2eSlice.begin(), mo2eSlice.end(), mo2e->begin());
      auto jop = make_shared<DimerJop>(nrasorb, system->norb(), make_shared<CSymMatrix>(nsysorb), mo2e);

      decimatetime.tick_print("make jop");

      map<BlockKey, shared_ptr<const Matrix>> hmap;
      map<BlockKey, shared_ptr<const Matrix>> spinmap;
      for (auto& c : civecs) {
        hmap.emplace(c.first, compute_sigma2e(c.second, jop));
        spinmap.emplace(c.first, compute_spin(c.second));
      }
      decimatetime.tick_print("compute renormalized 2e energy and spin");

      GammaForestProdASD forest(civecs);
      decimatetime.tick_print("construct GammaForestProdASD");

      forest.compute();
      decimatetime.tick_print("renormalize blocks");

      auto out = make_shared<DMRG_Block1>(move(forest), hmap, spinmap, ref->coeff()->slice_copy(ref->nclosed(), ref->nclosed()+ref->nact())->merge(system->coeff()));
      return out;
    }
  }
}

map<BlockKey, shared_ptr<const RASDvec>> RASD::diagonalize_site_RDM(const vector<shared_ptr<ProductRASCivec>>& civecs, const double perturbation) const {
  assert(civecs.size()==nstate_);

  Timer diagtime;

  // store the coefficients and sector bases by block: BlockKey corresponds to the block's information
  map<BlockKey, shared_ptr<Matrix>> coeffmap;

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
    // phase comes from rearranging |CI>|L> to |L>|CI>
    for (int ist = 0; ist < nstate_; ++ist) {
      for (auto& isec : civecs[ist]->sectors()) {
        const int nele_block = isec.first.nelea + isec.first.neleb;
        const int nele_ci = isec.second->det()->nelea() + isec.second->det()->neleb();
        if (detmap.find(isec.first)==detmap.end())
          detmap[isec.first] = isec.second->det();
        if ((nele_block*nele_ci)%2==1) {
          auto tmp = isec.second->copy();
          tmp->scale(-1.0);
          outer_products[isec.first].emplace_back(weights_[ist], tmp);
        }
        else {
          outer_products[isec.first].emplace_back(weights_[ist], isec.second);
        }
      }
    }

    diagtime.tick_print("organize outer products");
    // add in perturbative correction
    if (perturbation >= perturb_thresh_) {
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

    diagtime.tick_print("apply perturbation");

    for (auto& basis_sector : outer_products) {
      const int bsize = accumulate(basis_sector.second.begin(), basis_sector.second.end(), 0, [] (int x, tuple<double, shared_ptr<const Matrix>> m) { return x+get<1>(m)->mdim(); });
      auto tmp_mat = make_shared<Matrix>(get<1>(basis_sector.second.front())->ndim(), bsize);
      int current = 0;
      for (auto& ib : basis_sector.second) {
        copy_n(get<1>(ib)->data(), get<1>(ib)->size(), tmp_mat->element_ptr(0, current));
        current += get<1>(ib)->mdim();
      }
      for (int j = 0; j < bsize; ++j) {
        const double norm = blas::dot_product(tmp_mat->element_ptr(0,j), tmp_mat->ndim(), tmp_mat->element_ptr(0,j));
        if (norm > 1.0e-30)
          blas::scale_n(1.0/sqrt(norm), tmp_mat->element_ptr(0,j), tmp_mat->ndim());
      }
#ifdef HAVE_MPI_H
      tmp_mat->synchronize();
#endif
      cibasis.emplace(basis_sector.first, tmp_mat);
    }
    diagtime.tick_print("build basis out of outer products");
  }

  // RDM is block diagonal by sector
  for (auto& isec : cibasis) {
    shared_ptr<const Matrix> sectorbasis = isec.second;

    // build overlap matrix
    Matrix overlap(*sectorbasis % *sectorbasis);
#ifdef HAVE_MPI_H
    overlap.synchronize();
#endif
    diagtime.tick_print("build overlap");

    // build rdm
    Matrix rdm = *overlap.clone();
    for (auto& op : outer_products[isec.first]) {
      Matrix tmp(*get<1>(op) % *sectorbasis);
      dgemm_("T", "N", rdm.ndim(), rdm.mdim(), tmp.ndim(), get<0>(op), tmp.data(), tmp.ndim(), tmp.data(), tmp.ndim(), 1.0, rdm.data(), rdm.ndim());
    }
#ifdef HAVE_MPI_H
    rdm.synchronize();
#endif
    diagtime.tick_print("build rdm");

    Matrix orthonormalize(*overlap.tildex(1.0e-10));

    if (orthonormalize.mdim() > 0) {
      auto best_states = make_shared<Matrix>(orthonormalize % rdm * orthonormalize);
#ifdef HAVE_MPI_H
      best_states->synchronize();
#endif
      VectorB eigs(best_states->ndim());
      best_states->diagonalize(eigs);
      best_states = make_shared<Matrix>(orthonormalize * *best_states);
#ifdef HAVE_MPI_H
      best_states->synchronize();
      mpi__->broadcast(eigs.data(), eigs.size(), 0);
#endif
      coeffmap.emplace(isec.first, best_states);
      for (int i = 0; i < eigs.size(); ++i)
        singular_values.emplace(eigs(i), make_tuple(isec.first, i));
      diagtime.tick_print("diagonalize rdm and organize results");
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

    shared_ptr<const Matrix> coeff = coeffmap[bk];
    shared_ptr<const Matrix> sectorbasis = cibasis[bk];
    shared_ptr<const RASDeterminants> det = detmap[bk];

    auto tmp = make_shared<RASCivec>(det);
    dgemv_("N", det->size(), sectorbasis->mdim(), 1.0, sectorbasis->data(), sectorbasis->ndim(), coeff->element_ptr(0, position), 1, 1.0, tmp->data(), 1);
#ifdef HAVE_MPI_H
    tmp->synchronize();
#endif

    output_vectors[get<0>(i->second)].push_back(tmp);
    partial_trace += i->first;
  }
  const double total_trace = accumulate(singular_values.begin(), singular_values.end(), 0.0,
                                          [] (double x, pair<double, tuple<BlockKey, int>> s) { return x + s.first; } );
  cout << "  discarded weights: " << setw(12) << setprecision(8) << scientific <<  total_trace - partial_trace << fixed << endl;

  // Process into Dvecs: These vectors will become blocks so now BlockKey should be a descriptor of the block vector
  map<BlockKey, shared_ptr<const RASDvec>> out;
  cout << "  o Renormalized blocks have" << endl;
  for (auto& dvec : output_vectors) {
    auto tmp = make_shared<RASDvec>(dvec.second);
    cout << "    + " << tmp->ij() << " states with (" << tmp->det()->nelea() << ", " << tmp->det()->neleb() << ")" << endl;
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
    const int na = sdet->nelea() + dele.first;
    const int nb = sdet->neleb() + dele.second;
    if (na >= 0 && na <= sdet->norb() && nb >= 0 && nb <= sdet->norb()) {
      tdet = sdet->clone(na, nb);
      detmap[Tkey] = tdet;
    }
  }

  if (tdet && tdet->size()!=0) {
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
#ifdef HAVE_MPI_H
      tmp->synchronize();
#endif
      outer_products[Tkey].emplace_back(weight, tmp);
    }
  }
}

map<BlockKey, vector<shared_ptr<ProductRASCivec>>> RASD::diagonalize_site_and_block_RDM(const vector<shared_ptr<ProductRASCivec>>& civecs, const double perturbation) const {
  assert(civecs.size()==nstate_);

  Timer rdmtime(2);

  // for convenience since this will come up a lot in this function
  using ProdVec = vector<shared_ptr<ProductRASCivec>>;

  // first organize by environment block (right block in DMRG_Block2 object)
  shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(civecs.front()->left());
  assert(doubleblock);

  // store the coefficients and sector bases by block: BlockKey corresponds to the block's information
  map<BlockKey, shared_ptr<const Matrix>> coeffmap;

  // arrange all the 'singular values' to get the best ones
  multimap<double, tuple<BlockKey, const int>> singular_values;

  // arrange all the determinant objects for later use
  map<BlockKey, shared_ptr<const RASDeterminants>> detmap;

  // non-orthogonal basis in which to do the diagonalization
  map<BlockKey, ProdVec> cibasis;

  // list of all the vectors that get coupled together in the density matrix
  // rho = \sum_n \sum_{x,y in outer_products[n]} |x> w_n <y|
  map<BlockKey, vector<tuple<double, ProdVec>>> outer_products;

  // construct the non-orthogonal basis and the outer product vectors
  {
    // first, collect all of the vectors that belong in the state vectors
    for (int ist = 0; ist < nstate_; ++ist) {
      map<BlockKey, ProdVec> splitvec = civecs[ist]->split();
      for (auto& i : splitvec) {
        for (auto& isec : i.second.front()->sectors())
          if (detmap.find(isec.first)==detmap.end())
            detmap[isec.first] = isec.second->det();

        outer_products[i.first].emplace_back(weights_[ist], move(i.second));
      }
    }

    rdmtime.tick_print("outer products set up");;

    // add in perturbative correction
    if (perturbation >= perturb_thresh_) {
      map<BlockKey, vector<tuple<double, ProdVec>>> tmp_outerproducts = outer_products;
      for (auto& op : tmp_outerproducts) {
        for (auto& i : op.second) {
          // "diagonal" perturbation: sum_{i,j} [ (i^dagger j)_alpha (i^dagger_j)_beta ]
          apply_perturbation(get<1>(i), op.first, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, get<0>(i)*perturbation, outer_products);
          apply_perturbation(get<1>(i), op.first, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}, get<0>(i)*perturbation, outer_products);

          // spinflip perturbations
          apply_perturbation(get<1>(i), op.first, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}, get<0>(i)*perturbation, outer_products);
          apply_perturbation(get<1>(i), op.first, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, get<0>(i)*perturbation, outer_products);

          // ET/HT perturbations
          apply_perturbation(get<1>(i), op.first, {GammaSQ::AnnihilateAlpha}, get<0>(i)*perturbation, outer_products);
          apply_perturbation(get<1>(i), op.first, {GammaSQ::AnnihilateBeta}, get<0>(i)*perturbation, outer_products);
          apply_perturbation(get<1>(i), op.first, {GammaSQ::CreateAlpha}, get<0>(i)*perturbation, outer_products);
          apply_perturbation(get<1>(i), op.first, {GammaSQ::CreateBeta}, get<0>(i)*perturbation, outer_products);
        }
      }

      rdmtime.tick_print("perturbation applied");;
    }

#if 0
    // semi-normalize outer product vectors and pull the norm into the weight
    map<BlockKey, vector<tuple<double, ProdVec>>> tmp_outerproducts = outer_products;
    for (auto& op : tmp_outerproducts) {

    }
#endif
    for (auto& basis_sector : outer_products) {
      ProdVec tmpvec;
      for (auto& i : basis_sector.second)
        tmpvec.insert(tmpvec.end(), get<1>(i).begin(), get<1>(i).end());
      cibasis.emplace(basis_sector.first, move(tmpvec));
    }
    rdmtime.tick_print("basis built");;
  }

  // RDM is block diagonal by sector
  for (auto& isec : cibasis) {
    ProdVec sectorbasis = isec.second;

    // build overlap matrix
    Matrix overlap(sectorbasis.size(), sectorbasis.size());
    for (int i = 0; i < sectorbasis.size(); ++i) {
      for (int j = 0; j < i; ++j)
        overlap(i,j) = overlap(j,i) = sectorbasis[i]->dot_product(*sectorbasis[j]);
      overlap(i,i) = sectorbasis[i]->dot_product(*sectorbasis[i]);
    }
#ifdef HAVE_MPI_H
    overlap.synchronize();
#endif

    // build rdm
    Matrix rdm = *overlap.clone();
    for (auto& op : outer_products[isec.first]) {
      Matrix tmp(get<1>(op).size(), sectorbasis.size());
      for (int i = 0; i < tmp.mdim(); ++i)
        for (int j = 0; j < tmp.ndim(); ++j)
          tmp(j,i) = sectorbasis[i]->dot_product(*get<1>(op)[j]);
      dgemm_("T", "N", rdm.ndim(), rdm.mdim(), tmp.ndim(), get<0>(op), tmp.data(), tmp.ndim(), tmp.data(), tmp.ndim(), 1.0, rdm.data(), rdm.ndim());
    }
#ifdef HAVE_MPI_H
    rdm.synchronize();
#endif

    Matrix orthonormalize(*overlap.tildex(1.0e-11));

    if (orthonormalize.mdim() > 0) {
#ifdef HAVE_MPI_H
      orthonormalize.synchronize();
#endif
      rdmtime.tick_print("ortho built");

      auto best_states = make_shared<Matrix>(orthonormalize % rdm * orthonormalize);
#ifdef HAVE_MPI_H
      best_states->synchronize();
#endif
      VectorB eigs(best_states->ndim());
      best_states->diagonalize(eigs);
      best_states = make_shared<Matrix>(orthonormalize * *best_states);
#ifdef HAVE_MPI_H
      best_states->synchronize();
      mpi__->broadcast(eigs.data(), eigs.size(), 0);
#endif
      coeffmap.emplace(isec.first, best_states);
      for (int i = 0; i < eigs.size(); ++i)
        singular_values.emplace(eigs(i), make_tuple(isec.first, i));
    }
    rdmtime.tick_print("rdm block diagonalized");;
  }

  map<BlockKey, ProdVec> out;

  // Pick the top singular values
  int nvectors = 0;
  double partial_trace = 0.0;
  for (auto i = singular_values.rbegin(); i != singular_values.rend(); ++i, ++nvectors) {
    if (nvectors==ntrunc_) break;
    BlockKey bk = get<0>(i->second);
    const int position = get<1>(i->second);

    shared_ptr<const Matrix> coeff = coeffmap[bk];
    const ProdVec& sectorbasis = cibasis[bk];

    auto tmp = sectorbasis.front()->clone();
    for (int j = 0; j < sectorbasis.size(); ++j)
      tmp->ax_plus_y(coeff->element(j, position), *sectorbasis[j]);

#ifdef HAVE_MPI_H
    tmp->synchronize();
#endif
    out[BlockKey(tmp->nelea(), tmp->neleb())].push_back(tmp);
    partial_trace += i->first;
  }

  cout << "  o Renormalized blocks have" << endl;
  for (auto& o : out)
    cout << "    + " << o.second.size() << " states with (" << o.first.nelea << ", " << o.first.neleb << ")" << endl;

  const double total_trace = accumulate(singular_values.begin(), singular_values.end(), 0.0,
                                          [] (double x, pair<double, tuple<BlockKey, int>> s) { return x + s.first; } );
  cout << "  discarded weights: " << setw(12) << setprecision(8) << scientific <<  total_trace - partial_trace << fixed << endl;

  return out;
}

void RASD::apply_perturbation(const vector<shared_ptr<ProductRASCivec>>& ccvec, const BlockKey cckey, vector<GammaSQ> oplist,
                        const double weight, map<BlockKey, vector<tuple<double, vector<shared_ptr<ProductRASCivec>>>>>& outer_products) const
{
  pair<int, int> dele(0, 0);
  for (auto& op : oplist) {
    if (is_alpha(op))
      dele.first += (is_creation(op) ? 1 : -1);
    else
      dele.second += (is_creation(op) ? 1 : -1);
  }

  shared_ptr<const DMRG_Block> dmrgblock = ccvec.front()->left();
  shared_ptr<RASSpace> rasspace = ccvec.front()->space();

  const int tnelea = ccvec.front()->nelea() + dele.first;
  const int tneleb = ccvec.front()->neleb() + dele.second;

  const int blocknelea = cckey.nelea - dele.first;
  const int blockneleb = cckey.neleb - dele.second;

  const int norb = rasspace->norb();

  if (tnelea<0 || tneleb<0 || blocknelea<0 || blockneleb<0)
    return;

  for (auto& cc : ccvec) {
    auto out = make_shared<ProductRASCivec>(rasspace, dmrgblock, tnelea, tneleb);
    if (out->size() > 0) {
      ApplyOperator apply;
      for (auto& target_iter : out->sectors()) {
        const BlockKey target_key = target_iter.first;
        shared_ptr<RASBlockVectors> target_sector = target_iter.second;
        // perturbation is applied only to the CI part
        if (cc->contains_block(target_key)) {
          shared_ptr<const RASBlockVectors> source_sector = cc->sector(target_key);
          assert(target_sector->mdim()==source_sector->mdim());
          if (oplist.size()==1) {
            for (int p = 0; p < norb; ++p)
              apply(1.0, *source_sector, *target_sector, oplist, {p});
          } else if (oplist.size()==2) {
            for (int p = 0; p < norb; ++p)
              for (int q = 0; q < norb; ++q)
                apply(1.0, *source_sector, *target_sector, oplist, {p, q});
          } else {
            assert(false);
          }
        }
      }
#ifdef HAVE_MPI_H
      out->synchronize();
#endif
      const BlockKey target_key(blocknelea, blockneleb);
      outer_products[target_key].emplace_back(weight, vector<shared_ptr<ProductRASCivec>>{{out}});
    }
  }
}
