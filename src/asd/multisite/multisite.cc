//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
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

#include <src/asd/multisite/multisite.h>
#include <src/scf/hf/fock.h>
#include <src/mat1e/overlap.h>
#include <src/wfn/localization.h>

using namespace std;
using namespace bagel;

// constructor
MultiSite::MultiSite(shared_ptr<const PTree> input, shared_ptr<const Reference> ref) : hf_ref_(ref) {
  cout << string(60, '=') << endl;
  cout << string(15, ' ') << "Constructing MultiSite" << endl;
  cout << string(60, '=') << endl;

  // localize closed and virtual orbitals (optional)
  shared_ptr<const PTree> localization_data = input->get_child_optional("localization");
  if (localization_data)
    localize(localization_data);
  else
    mref_ = make_shared<Reference>(*hf_ref_);
}


void MultiSite::localize(shared_ptr<const PTree> localization_data) {

  if (!localization_data->get_child_optional("region_sizes")) throw runtime_error("region_sizes has to be provided to do localization");
  vector<int> region_sizes = localization_data->get_vector<int>("region_sizes");
  assert(accumulate(region_sizes.begin(), region_sizes.end(), 0) == hf_ref_->geom()->natom());

  string localizemethod = localization_data->get<string>("algorithm", "pm");

  shared_ptr<OrbitalLocalization> localization;
  auto input_data = make_shared<PTree>(*localization_data);
  input_data->erase("virtual"); input_data->put("virtual", true);

  if (localizemethod == "region") {
    localization = make_shared<RegionLocalization>(input_data, hf_ref_);
  } else if (localizemethod == "pm" || localizemethod == "pipek-mezey") {
    input_data->erase("type"); input_data->put("type", "region");
    localization = make_shared<PMLocalization>(input_data, hf_ref_);
  } else throw runtime_error("Unrecognized orbital localization method");

  shared_ptr<const Matrix> local_coeff = localization->localize();
  vector<pair<int, int>> orbital_subspaces = localization->orbital_subspaces();
  assert(orbital_subspaces.size() == 2); // closed and active only

  vector<pair<int, int>> region_bounds;
  {
    int atom_start = 0;
    int current = 0;
    for (int natom : region_sizes) {
      const int bound_start = current;
      for (int iatom = 0; iatom != natom; ++iatom)
        current += hf_ref_->geom()->atoms(atom_start+iatom)->nbasis();
      region_bounds.emplace_back(bound_start, current);
      atom_start += natom;
    }
  }

  assert(hf_ref_->coeff()->mdim() == accumulate(orbital_subspaces.begin(), orbital_subspaces.end(), 0ull,
                               [] (int o, const pair<int, int>& p) { return o + p.second - p.first; }));

  Matrix ShalfC = static_cast<Matrix>(Overlap(hf_ref_->geom()));
  ShalfC.sqrt();
  ShalfC *= *local_coeff;

  shared_ptr<const Fock<1>> fock;
  {
    shared_ptr<const Matrix> density = hf_ref_->coeff()->form_density_rhf(hf_ref_->nclosed());
    const MatView ccoeff = hf_ref_->coeff()->slice(0, hf_ref_->nclosed());
    fock = make_shared<const Fock<1>>(hf_ref_->geom(), hf_ref_->hcore(), density, ccoeff);
  }

  const int nsites = region_sizes.size();
  auto out_coeff = local_coeff->clone();

  for (const pair<int, int>& subspace : orbital_subspaces) {
    vector<set<int>> orbital_sets(nsites);
    const int nsuborbs = subspace.second - subspace.first;

    vector<vector<double>> lowdin_populations(nsuborbs);
    for (auto& p : lowdin_populations)
      p.resize(nsites);

    Matrix Q(nsuborbs, nsuborbs);
    for (int site = 0; site != nsites; ++site) {
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
      imo += subset.size();
    }
  }

  mref_ = make_shared<Reference>(*hf_ref_, make_shared<Coeff>(move(*out_coeff)));
}


