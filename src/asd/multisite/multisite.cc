//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
// Copyright (C) 2014 Shane Parker
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

#include <src/asd/multisite/multisite.h>
#include <src/scf/hf/fock.h>
#include <src/mat1e/overlap.h>
#include <src/wfn/localization.h>
#include <src/util/io/moldenout.h>

using namespace std;
using namespace bagel;

// constructor
MultiSite::MultiSite(shared_ptr<const PTree> input, shared_ptr<const Reference> ref, const int nsites) : input_(input), hf_ref_(ref), nsites_(nsites) {
  cout << string(60, '=') << endl;
  cout << string(15, ' ') << "Construct MultiSite" << endl;
  cout << string(60, '=') << endl;

  charge_ = input_->get<int>("charge", 0);
  nspin_ = input_->get<int>("nspin", 0);
  active_electrons_ = input_->get_vector<int>("active_electrons");
  active_sizes_ = input_->get_vector<int>("active_sizes");
  region_sizes_ = input_->get_vector<int>("region_sizes");
  assert(accumulate(region_sizes_.begin(), region_sizes_.end(), 0) == hf_ref_->geom()->natom());
}


void MultiSite::compute() {
  // construct Fock matrix for localization and canonicalization
  shared_ptr<const Fock<1>> fock;
  {
    shared_ptr<const Matrix> density = hf_ref_->coeff()->form_density_rhf(hf_ref_->nclosed());
    const MatView ccoeff = hf_ref_->coeff()->slice(0, hf_ref_->nclosed());
    fock = make_shared<const Fock<1>>(hf_ref_->geom(), hf_ref_->hcore(), density, ccoeff);
  }

  // localize closed and virtual orbitals (optional)
  shared_ptr<const PTree> localization_data = input_->get_child_optional("localization");
  if (localization_data) {
    localize(localization_data, fock);
    MoldenOut mfile("localized.molden");
    mfile << sref_->geom();
    mfile << sref_;
  } else {
    sref_ = make_shared<Reference>(*hf_ref_);
  }

  // two options : 1) manually assign active orbitals to subspaces  2) use projection to do it automatically
  set_active_orbitals();

  // canonicalize active orbitals within each subspace
  canonicalize(fock);

}


void MultiSite::localize(shared_ptr<const PTree> localize_data, shared_ptr<const Matrix> fock) {
  string localizemethod = localize_data->get<string>("algorithm", "pm");

  shared_ptr<OrbitalLocalization> localization;
  auto input_data = make_shared<PTree>(*localize_data);
  input_data->erase("virtual"); input_data->put("virtual", true);
  if (input_data->get_child_optional("region_sizes")) {
    cout << "WARNING : The region_sizes keyword in localization input will be overwritten by that from MultiSite." << endl;
    input_data->erase("region_sizes");
  }

  if (localizemethod == "region") {
    localization = make_shared<RegionLocalization>(input_data, hf_ref_, region_sizes_);
  } else if (localizemethod == "pm" || localizemethod == "pipek-mezey") {
    input_data->erase("type"); input_data->put("type", "region");
    localization = make_shared<PMLocalization>(input_data, hf_ref_, region_sizes_);
  } else throw runtime_error("Unrecognized orbital localization method");

  shared_ptr<const Matrix> local_coeff = localization->localize();
  vector<pair<int, int>> orbital_subspaces = localization->orbital_subspaces();
  assert(orbital_subspaces.size() == 2); // closed and active only

  vector<pair<int, int>> region_bounds;
  {
    int atom_start = 0;
    int current = 0;
    for (int natom : region_sizes_) {
      const int bound_start = current;
      for (int iatom = 0; iatom != natom; ++iatom)
        current += hf_ref_->geom()->atoms(atom_start+iatom)->nbasis();
      region_bounds.emplace_back(bound_start, current);
      atom_start += natom;
    }
  }
  assert(region_bounds.size() == nsites_);

  vector<vector<pair<int, int>>> site_bounds;

  assert(hf_ref_->coeff()->mdim() == accumulate(orbital_subspaces.begin(), orbital_subspaces.end(), 0ull,
                               [] (int o, const pair<int, int>& p) { return o + p.second - p.first; }));

  Matrix ShalfC = static_cast<Matrix>(Overlap(hf_ref_->geom()));
  ShalfC.sqrt();
  ShalfC *= *local_coeff;

  auto out_coeff = local_coeff->clone();
  
  for (const pair<int, int>& subspace : orbital_subspaces) {
    vector<set<int>> orbital_sets(nsites_);
    const int nsuborbs = subspace.second - subspace.first;

    vector<pair<int, int>> subspace_site_bounds;

    vector<vector<double>> lowdin_populations(nsuborbs);
    for (auto& p : lowdin_populations)
      p.resize(nsites_);

    Matrix Q(nsuborbs, nsuborbs);
    for (int site = 0; site != nsites_; ++site) {
      const pair<int, int> bounds = region_bounds[site];
      const int nregbasis = bounds.second - bounds.first;
      dgemm_("T", "N", nsuborbs, nsuborbs, nregbasis, 1.0, ShalfC.element_ptr(bounds.first, subspace.first), ShalfC.ndim(),
                                                           ShalfC.element_ptr(bounds.first, subspace.first), ShalfC.ndim(), 
                                                      0.0, Q.data(), Q.ndim());
      for (int orb = 0; orb != nsuborbs; ++orb)
        lowdin_populations[orb][site] = Q(orb, orb) * Q(orb, orb);
    }

    for (int orb = 0; orb != nsuborbs; ++orb) {
      //pick the largest value to assign it to sites
      vector<double>& pops = lowdin_populations[orb];
      auto maxiter = max_element(pops.begin(), pops.end());
      const int maxsite = maxiter - pops.begin();
      orbital_sets[maxsite].insert(orb + subspace.first);
    }

    size_t imo = subspace.first;

    for (auto& subset : orbital_sets) {
      if (subset.empty()) 
        throw runtime_error("one of the active subspaces has bad definition, no matching active orbitals are localized on it, please redefine subspaces");
      Matrix subspace(out_coeff->ndim(), subset.size());
      int pos = 0;
      for (const int& i : subset)
        copy_n(local_coeff->element_ptr(0, i), local_coeff->ndim(), subspace.element_ptr(0, pos++));

      Matrix subfock(subspace % *fock * subspace);
      VectorB eigs(subspace.ndim());
      subfock.diagonalize(eigs);
      subspace *= subfock;

      copy_n(subspace.data(), subspace.size(), out_coeff->element_ptr(0, imo));
      subspace_site_bounds.emplace_back(imo, imo+subset.size());
      imo += subset.size();
    }
    site_bounds.emplace_back(move(subspace_site_bounds));
  }

  assert(site_bounds.front().front().first < site_bounds.back().front().first);
  
  sref_ = make_shared<Reference>(hf_ref_->geom(), make_shared<Coeff>(move(*out_coeff)), hf_ref_->nclosed(), hf_ref_->nact(), hf_ref_->nvirt());
}


void MultiSite::set_active_orbitals() {
  
  // collect orbital subspaces info
  const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
  const int nclosed = (hf_ref_->geom()->nele() - charge_ - nactele) / 2;
  assert((hf_ref_->geom()->nele() - charge_ - nactele) % 2 == 0);
  const int nactive = accumulate(active_sizes_.begin(), active_sizes_.end(), 0);
  const int nvirt = hf_ref_->coeff()->mdim() - nclosed - nactive;
  const int multisitebasis = hf_ref_->geom()->nbasis();

  auto active_subspaces = input_->get_child("active_subspaces");
  set<int> active_set;

  auto out_coeff = sref_->coeff()->clone();

  if (active_subspaces->size() == nsites_) {
    // active subspaces info is provided, just reorder the orbitals
    
    vector<set<int>> active_orbitals;
    for (auto site : *active_subspaces) {
      vector<int> active_orbs = site->get_vector<int>("");
      for_each(active_orbs.begin(), active_orbs.end(), [](int& x) { --x; });
      active_orbitals.emplace_back(active_orbs.begin(), active_orbs.end());
      active_set.insert(active_orbs.begin(), active_orbs.end());
    }
    assert(active_set.size() == nactive);
    
    auto in_coeff = sref_->coeff();
    assert(in_coeff->ndim() == multisitebasis);
  
    int active_position = nclosed;
    for (auto site : active_orbitals) {
      for (int aorb : site)
        copy_n(in_coeff->element_ptr(0, aorb), multisitebasis, out_coeff->element_ptr(0, active_position++));
    }
    assert(active_position == nclosed + nactive);
    
  } else if (active_subspaces->size() == 1){
    // if no manually assigned active orbital subspaces, do projection

    for (auto site : *active_subspaces) {
      vector<int> active_vec = site->get_vector<int>("");
      for_each(active_vec.begin(), active_vec.end(), [](int& x) { --x; });
      active_set.insert(active_vec.begin(), active_vec.end());
    }
    assert(active_set.size() == nactive);

    // delocalized active orbital set
    auto act_orbs = make_shared<Matrix>(multisitebasis, nactive);
    int pos = 0;
    for (auto iorb : active_set)
      copy_n(sref_->coeff()->element_ptr(0, iorb), multisitebasis, act_orbs->element_ptr(0, pos++));
  
    // region bound info
    vector<pair<int, int>> region_bounds;
    int atomstart = 0; int basis_start = 0;
    for (int isize : region_sizes_) {
      int nbasis = 0;
      for (int iatom = atomstart; iatom != atomstart+isize; ++iatom)
        nbasis += sref_->geom()->atoms(iatom)->nbasis();
      region_bounds.emplace_back(basis_start, basis_start+nbasis);
      atomstart += isize;
      basis_start += nbasis;
    }
  
    // do SVD and find transformation matrix
    Matrix transform(nactive, nactive);
    {
      int ipos = 0;
      for (int isite = 0; isite != nsites_; ++isite) {
        pair<int, int> bounds = region_bounds[isite];
        const int ntot = nactive;
        const int nsub = active_sizes_[isite];
        auto actorb = act_orbs->get_submatrix(bounds.first, 0, bounds.second-bounds.first, ntot);
        Matrix SVD(*actorb % *actorb);
        VectorB eigs(SVD.ndim());
        SVD.diagonalize(eigs);
        const MatView sub_trans = SVD.slice(ntot-nsub, ntot);
        copy_n(sub_trans.data(), sub_trans.size(), transform.element_ptr(0, ipos));
        ipos += nsub;
      }
      assert(ipos == nactive);
      // Lowdin orthonormalization
      Matrix inv_hf(transform % transform);
      inv_hf.inverse_half();
      transform *= inv_hf;
      cout << endl << string(6, '*') << "  If linear dependency is detected, you have to find better active orbitals to do projection  " << string(6, '*') << endl << endl;
    }
    
    // project coeff
    Matrix projected(*act_orbs * transform);
    copy_n(projected.data(), projected.size(), out_coeff->element_ptr(0, nclosed));

  } else {
    throw runtime_error("Must either provide one total active orbital list or assign active orbitals for each fragment");
  }

  int closed_position = 0;
  int virt_position = nclosed + nactive;
  for (int iorb = 0; iorb != out_coeff->mdim(); ++iorb)
    if (active_set.count(iorb) == 0)
      copy_n(sref_->coeff()->element_ptr(0, iorb), multisitebasis, out_coeff->element_ptr(0, (closed_position < nclosed ? closed_position++ : virt_position++)));
  assert(closed_position == nclosed);
  assert(virt_position == out_coeff->mdim());
  
  sref_ = make_shared<Reference>(sref_->geom(), make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, nvirt);
}


void MultiSite::canonicalize(shared_ptr<const Matrix> fock) {
  //canonicalize active orbitals within subspace
  const int nclosed = sref_->nclosed();
  
  // rotate within each active subspace
  const MatView active_mos = sref_->coeff()->slice(nclosed, nclosed+sref_->nact());
  Matrix fock_mo(active_mos % *fock * active_mos);
  VectorB eigs(active_mos.mdim());
  shared_ptr<const Matrix> active_transformation = fock_mo.diagonalize_blocks(eigs, active_sizes_);
  const Matrix transformed_mos(active_mos * *active_transformation);
  shared_ptr<Matrix> scoeff = sref_->coeff()->copy();
  copy_n(transformed_mos.data(), transformed_mos.size(), scoeff->element_ptr(0, nclosed));

  sref_ = make_shared<Reference>(*sref_, make_shared<Coeff>(move(*scoeff)));
}


shared_ptr<Reference> MultiSite::build_reference(const int site, const vector<bool> meanfield) const {
  assert(meanfield.size() == nsites_ && site < nsites_ && site >=0);

  vector<shared_ptr<const MatView>> closed_orbitals = { make_shared<MatView>(sref_->coeff()->slice(0, sref_->nclosed())) };
  const int act_start = accumulate(active_sizes_.begin(), active_sizes_.begin()+site, sref_->nclosed());
  const int act_fence = act_start + active_sizes_.at(site);
  const MatView active_orbitals = sref_->coeff()->slice(act_start, act_fence);

  int current = sref_->nclosed();
  for (int i = 0; i != nsites_; ++i) {
    if (meanfield[i] && i != site)
      closed_orbitals.push_back(make_shared<const MatView>(sref_->coeff()->slice(current, current+(active_electrons_.at(i)+1)/2)));
    current += active_sizes_.at(i);
  }

  const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [](int x, shared_ptr<const MatView>m) { return x + m->mdim(); });
  const int nact = active_orbitals.mdim();

  auto out = make_shared<Matrix>(sref_->geom()->nbasis(), nclosed+nact);

  current = 0;
  closed_orbitals.push_back(make_shared<MatView>(active_orbitals));
  for (auto& orbitals : closed_orbitals) {
    copy_n(orbitals->data(), orbitals->size(), out->element_ptr(0, current));
    current += orbitals->mdim();
  }

  return make_shared<Reference>(sref_->geom(), make_shared<Coeff>(move(*out)), nclosed, nact, 0);
}


