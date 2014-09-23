//
// BAGEL - Parallel electron correlation program.
// Filename: multisite_scf.cc
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

#include <src/multisite/multisite.h>
#include <src/scf/fock.h>
#include <src/wfn/construct_method.h>
#include <src/scf/scf.h>

#include <src/molecule/localization.h>
#include <src/molecule/overlap.h>

using namespace std;
using namespace bagel;

/// Localization of orbitals onto site fragments --
///  EVERYTHING (including virtuals) is localized.
///  In principle, ambiguous orbitals could be handled but they are not yet.
void MultiSite::localize(const shared_ptr<const PTree> idata, shared_ptr<const Matrix> fock) {
  string localizemethod = idata->get<string>("algorithm", "pm");

  shared_ptr<OrbitalLocalization> localization;
  auto input_data = make_shared<PTree>(*idata);
  input_data->erase("virtual"); input_data->put("virtual", "true");
  if (input_data->get_child_optional("region_sizes")) {
    cout << "WARNING: The region_sizes keyword will be ignored for multisite preparation." << endl;
    input_data->erase("region_sizes");
  }

  vector<int> sizes;
  for (auto& g : geoms_)
    sizes.push_back(g->natom());

  if (localizemethod == "region") {
    localization = make_shared<RegionLocalization>(input_data, sref_, sizes);
  }
  else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey") {
    input_data->erase("type"); input_data->put("type", "region");
    localization = make_shared<PMLocalization>(input_data, sref_, sizes);
  }
  else throw std::runtime_error("Unrecognized orbital localization method");

  shared_ptr<const Matrix> local_coeff = localization->localize();
  vector<pair<int, int>> orbital_subspaces = localization->orbital_subspaces();

  // this function requires exactly two subspaces (closed and virtual) to have been localized
  assert(orbital_subspaces.size()==2);

  vector<pair<int, int>> region_bounds;
  {
    int current = 0;
    for (auto& g : geoms_) {
      region_bounds.emplace_back(current, current + g->nbasis());
      current += g->nbasis();
    }
  }

  vector<vector<pair<int, int>>> site_bounds;

  assert(sref_->coeff()->mdim() == accumulate(orbital_subspaces.begin(), orbital_subspaces.end(), 0ull,
                             [] (unsigned long long o, const pair<int,int>& p) { return o + p.second - p.first; }));

  shared_ptr<const Geometry> sgeom = sref_->geom();

  Matrix ShalfC = static_cast<Matrix>(Overlap(sgeom));
  ShalfC.sqrt();
  ShalfC *= *local_coeff;

  auto out_coeff = local_coeff->copy();

  for (const pair<int, int>& subspace : orbital_subspaces) {
    vector<set<int>> orbital_sets(nsites_);
    const int nsuborbs = subspace.second - subspace.first;

    vector<pair<int,int>> subspace_site_bounds;

    vector<vector<double>> lowdin_populations(nsuborbs);
    for (auto& p : lowdin_populations)
      p.resize(nsites_);

    Matrix Q(nsuborbs, nsuborbs);
    for (int site = 0; site < nsites_; ++site) {
      const pair<int, int> bounds = region_bounds[site];
      const int nregbasis = bounds.second - bounds.first;
      dgemm_("T", "N", nsuborbs, nsuborbs, nregbasis, 1.0, ShalfC.element_ptr(bounds.first, subspace.first), ShalfC.ndim(),
                                                           ShalfC.element_ptr(bounds.first, subspace.first), ShalfC.ndim(),
                                                      0.0, Q.data(), Q.ndim());

      for (int orb = 0; orb < nsuborbs; ++orb)
        lowdin_populations[orb][site] = Q(orb,orb) * Q(orb,orb);
    }

    for (int orb = 0; orb < nsuborbs; ++orb) {
      // pick the largest value to assign it to sites
      vector<double>& pops = lowdin_populations[orb];
      auto maxiter = max_element(pops.begin(), pops.end());
      const int maxsite = maxiter - pops.begin();
      orbital_sets[maxsite].insert(orb + subspace.first);
    }

    size_t imo = subspace.first;

    for (auto& subset : orbital_sets) {
      if (subset.empty()) continue;
      Matrix subspace(out_coeff->ndim(), subset.size());
      int pos = 0;
      for (const int& i : subset)
        copy_n(local_coeff->element_ptr(0, i), local_coeff->ndim(), subspace.element_ptr(0, pos++));

      Matrix subfock(subspace % *fock * subspace);
      VectorB eigs(subspace.mdim());
      subfock.diagonalize(eigs);
      subspace = subspace * subfock;

      copy_n(subspace.data(), subspace.size(), out_coeff->element_ptr(0, imo));
      subspace_site_bounds.emplace_back(imo, imo+subset.size());
      imo += subset.size();
    }
    site_bounds.emplace_back(move(subspace_site_bounds));
  }

  // Beware: if the number of subspaces changes beyond just closed and virtual, this will need to be updated
  if (site_bounds.front().front().first < site_bounds.back().front().first) {
    closed_bounds_ = site_bounds.front();
    virt_bounds_ = site_bounds.back();
  }
  else {
    closed_bounds_ = site_bounds.back();
    virt_bounds_ = site_bounds.front();
  }
  active_bounds_.clear();

  sref_ = make_shared<Reference>(*sref_, make_shared<Coeff>(move(*out_coeff)));
}

///  Separate SVDs are done for each monomer
void MultiSite::set_active(const std::shared_ptr<const PTree> idata) {
  // TODO needs clean up
  auto active = idata->get_child_optional("active");
  if (!active)
    throw logic_error("\"active\" keyword must be set for multi-site algorithm.");

  if (active->size()!=1 && active->size()!=nsites_)
    throw logic_error("Either one or all of the sites must have active spaces specified");

  vector<set<int>> active_orbitals;
  for (auto site : *active) {
    vector<int> active_orbs = site->get_vector<int>("");
    for_each(active_orbs.begin(), active_orbs.end(), [] (int& x) { x--; }); // subtract one to make indexing match
    active_orbitals.emplace_back(active_orbs.begin(), active_orbs.end());
  }
  active_orbitals.resize(nsites_, active_orbitals.front());

  // Make new References
  for (int site = 0; site < nsites_; ++site)
    active_refs_.push_back(isolated_refs_[site]->set_active(active_orbitals[site]));

  const int nclosed = accumulate(active_refs_.begin(), active_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nclosed(); });
  const int nactive = accumulate(active_refs_.begin(), active_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nact(); });
  const int nvirt = accumulate(active_refs_.begin(), active_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nvirt(); });

  // TODO: this implementation requires specifying the number of active orbitals that are coming from each subset.
  //  This is probably fine, but it is not strictly necessary.

  // tuple:
  //  matrix --> reference active orbitals
  //  pair   --> bounds on subspace that can be mixed to produce localized orbitals
  //  int    --> number of orbitals to form
  //  string --> name of subspace (just for pretty printing)
  //  bool   --> closed(true)/virtual(false)
  vector<tuple<shared_ptr<const MatView>, pair<int, int>, int, string, bool>> svd_info;

  for (int site = 0; site < nsites_; ++site) {
    shared_ptr<const Reference> isoref = isolated_refs_[site];
    shared_ptr<const Reference> actref = active_refs_[site];
    const int nclosed_occ = isoref->nclosed() - actref->nclosed();
    const int nvirt_occ = isoref->nvirt() - actref->nvirt();
    auto actview = make_shared<const MatView>(actref->coeff()->slice(actref->nclosed(), actref->nclosed()+actref->nact()));
    svd_info.emplace_back(actview, closed_bounds_[site], nclosed_occ, to_string(site), true);
    svd_info.emplace_back(actview, virt_bounds_[site], nvirt_occ, to_string(site), false);
  }

  Overlap S(sref_->geom());

  shared_ptr<Matrix> out_coeff = sref_->coeff()->copy();
  size_t closed_position = 0;
  size_t active_position = nclosed;
  size_t virt_position = nclosed + nactive;

  for (auto& subset : svd_info) {
    const MatView& active = *get<0>(subset);
    pair<int, int> bounds = get<1>(subset);
    const int norb = get<2>(subset);
    const string site_name = get<3>(subset);
    const bool closed = get<4>(subset);

    const MatView subcoeff = sref_->coeff()->slice(bounds.first, bounds.second);

    const Matrix overlaps(active % S * subcoeff);

    multimap<double, int> norms;

    for(int i = 0; i < overlaps.mdim(); ++i) {
      const double norm = blas::dot_product(overlaps.element_ptr(0, i), overlaps.ndim(), overlaps.element_ptr(0, i));
      norms.emplace(norm, i);
    }

    active_thresh_ = input_->get<double>("active_thresh", 0.5);
    cout << endl << "  o Forming dimer's active orbitals arising from " << (closed ? "closed " : "virtual ") << "orbitals of site "
                 << site_name << ". Threshold for inclusion in cadidate space: " << setw(6) << setprecision(3) << active_thresh_ << endl;

    vector<int> active_list;
    double max_overlap, min_overlap;
    {
      auto end = norms.rbegin(); advance(end, norb);
      end = find_if(end, norms.rend(), [this] (const pair<const double, int>& p) { return p.first < active_thresh_; });
      for_each(norms.rbegin(), end, [&active_list] (const pair<const double, int>& p) { active_list.emplace_back(p.second); });
      auto mnmx = minmax_element(norms.rbegin(), end);
      tie(min_overlap, max_overlap) = make_tuple(mnmx.first->first, mnmx.second->first);
    }

    const int active_size = active_list.size();
    cout << "    - size of candidate space: " << active_size << endl;
    cout << "    - largest overlap with monomer space: " << max_overlap << ", smallest: " << min_overlap << endl;

    if (active_size != norb) {
      cout << "  o Performing SVD in candidate space" << endl;
      Matrix subspace(out_coeff->ndim(), active_size);

      int ii = 0;
      for (int& i : active_list)
        copy_n(subcoeff.element_ptr(0, i), subcoeff.ndim(), subspace.element_ptr(0, ii++));

      Matrix Sactive(active % S * active);
      Sactive.inverse_half();

      Matrix projector( Sactive * ( active % S * subspace ) );
      VectorB singulars(active_size);
      shared_ptr<Matrix> Vt;
      tie(ignore, Vt) = projector.svd(singulars.data());

      cout << "    - largest singular value: " << singulars[0] << ", smallest: " << singulars[norb-1] << endl;
      cout << "    - norb: " << norb << ", sum of highest singular values: " << accumulate(singulars.begin(), singulars.begin()+norb, 0.0) << endl;

      subspace = subspace ^ *Vt;

      for (size_t i = 0; i < subcoeff.mdim(); ++i) {
        if (count(active_list.begin(),active_list.end(),i)==0)
          copy_n(subcoeff.element_ptr(0, i), subcoeff.ndim(), out_coeff->element_ptr(0, ( closed ? closed_position++ : virt_position++ )));
      }
      copy_n(subspace.data(), subspace.size(), out_coeff->element_ptr(0, active_position));
      active_position += norb;

      for (size_t i = norb; i < active_size; ++i)
        copy_n(subspace.element_ptr(0, i), subspace.ndim(), out_coeff->element_ptr(0, ( closed ? closed_position++ : virt_position++ )));
    }
    else {
      set<int> active_set(active_list.begin(), active_list.end());
      for (size_t i = 0; i < subcoeff.mdim(); ++i)
        if (active_set.count(i)==0)
          copy_n(subcoeff.element_ptr(0, i), subcoeff.ndim(), out_coeff->element_ptr(0, ( closed ? closed_position++ : virt_position++ )));
        else
          copy_n(subcoeff.element_ptr(0, i), subcoeff.ndim(), out_coeff->element_ptr(0, active_position++));
    }
  }

  sref_ = make_shared<Reference>(sref_->geom(), make_shared<Coeff>(*out_coeff), nclosed, nactive, nvirt);
}

// RHF and then localize
void MultiSite::scf(const shared_ptr<const PTree> idata) {
  Timer stopwatch;

  // SCF
  auto hfdata = idata->get_child_optional("hf") ? idata->get_child_optional("hf") : make_shared<PTree>();
  shared_ptr<SCF> rhf = dynamic_pointer_cast<SCF>(construct_method("hf", hfdata, sref_->geom(), sref_));
  rhf->compute();
  sref_ = rhf->conv_to_ref();
  stopwatch.tick_print("Multisite SCF");

  const int nclosed = sref_->nclosed();

  shared_ptr<const Matrix> density = sref_->coeff()->form_density_rhf(nclosed);
  const MatView coeff = sref_->coeff()->slice(0,nclosed);

  shared_ptr<const PTree> localize_data = idata->get_child_optional("localization");
  if (!localize_data) localize_data = make_shared<const PTree>();

  auto fock  = make_shared<const Fock<1>>(sref_->geom(), sref_->hcore(), density, coeff);
  stopwatch.tick_print("Multisite Fock matrix formation");

  localize(localize_data, fock);
  stopwatch.tick_print("Site localization");

  set_active(idata);
  stopwatch.tick_print("Set active space");

  MatView active_mos = sref_->coeff()->slice(sref_->nclosed(), sref_->nclosed() + sref_->nact());
  Matrix fock_mo(active_mos % *fock * active_mos);
  VectorB eigs(active_mos.mdim());
  vector<int> active_blocks;
  for (auto& r : active_refs_)
    active_blocks.push_back(r->nact());
  shared_ptr<Matrix> active_transformation = fock_mo.diagonalize_blocks(eigs, active_blocks);
  Matrix transformed_mos(active_mos * *active_transformation);
  shared_ptr<Matrix> scoeff = sref_->coeff()->copy();
  scoeff->copy_block(0, sref_->nclosed(), scoeff->ndim(), transformed_mos.mdim(), transformed_mos);
  sref_ = make_shared<Reference>(*sref_, make_shared<Coeff>(move(*scoeff)));
}
