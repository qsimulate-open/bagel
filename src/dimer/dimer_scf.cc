//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_scf.cc
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

#include <tuple>

#include <src/dimer/dimer.h>
#include <src/scf/fock.h>
#include <src/molecule/localization.h>
#include <src/util/string.h>

using namespace std;
using namespace bagel;

/// Localization of dimer orbitals --
///  Prelocalize flag determines whether localization is before or after assigning orbitals to active space.
///  The main difference is that EVERYTHING (including virtuals) is localized for prelocalize and it's okay
///  in principle for ambiguous orbitals in the prelocalize (they just can't become active orbitals).
void Dimer::localize(const shared_ptr<const PTree> idata, shared_ptr<const Matrix> fock, const bool localize_first) {
  string localizemethod = idata->get<string>("algorithm", "pm");

  shared_ptr<OrbitalLocalization> localization;
  auto input_data = make_shared<PTree>(*idata);
  if (localize_first) { input_data->erase("virtual"); input_data->put("virtual", "true"); }
  vector<int> sizes = { geoms_.first->natom(), geoms_.second->natom() };
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

  vector<set<int>> subsets_A;
  vector<set<int>> subsets_B;
  vector<set<int>> ambiguous_subsets;

  pair<int, int> nbasis = make_pair(geoms_.first->nbasis(), geoms_.second->nbasis());

  if (localize_first)
    assert(sref_->coeff()->mdim() == accumulate(orbital_subspaces.begin(), orbital_subspaces.end(), 0ull,
                               [] (unsigned long long o, const pair<int,int>& p) { return o + p.second - p.first; }));

  Matrix ShalfC = static_cast<Matrix>(Overlap(sgeom_));
  ShalfC.sqrt();
  ShalfC *= *local_coeff;

  for (const pair<int, int>& subspace : orbital_subspaces) {
    set<int> set_A, set_B, ambiguous;
    const int nsuborbs = subspace.second - subspace.first;

    vector<double> lowdin_A(nsuborbs, 0.0);
    vector<double> lowdin_B(nsuborbs, 0.0);

    {
      Matrix Q_A(nsuborbs, nsuborbs);
      dgemm_("T", "N", nsuborbs, nsuborbs, nbasis.first, 1.0, ShalfC.element_ptr(0,subspace.first), ShalfC.ndim(),
                                                              ShalfC.element_ptr(0,subspace.first), ShalfC.ndim(),
                                                         0.0, Q_A.data(), Q_A.ndim());
      Matrix Q_B(nsuborbs, nsuborbs);
      dgemm_("T", "N", nsuborbs, nsuborbs, nbasis.second, 1.0, ShalfC.element_ptr(nbasis.first,subspace.first), ShalfC.ndim(),
                                                               ShalfC.element_ptr(nbasis.first,subspace.first), ShalfC.ndim(),
                                                          0.0, Q_B.data(), Q_B.ndim());
      for (int i = 0; i < nsuborbs; ++i) {
        lowdin_A[i] = Q_A(i,i) * Q_A(i,i);
        lowdin_B[i] = Q_B(i,i) * Q_B(i,i);
      }
    }

    for (int i = subspace.first; i < subspace.second; ++i) {
      const double Anorm = lowdin_A[i - subspace.first];
      const double Bnorm = lowdin_B[i - subspace.first];
      const double polarized = (Anorm - Bnorm)/(Anorm + Bnorm);

      if (polarized > 0.5)
        set_A.insert(i);
      else if (polarized < -0.5)
        set_B.insert(i);
      else {
        ambiguous.insert(i);
        // in principle, ambiguous orbitals could be handled if localizing first. Not implemented yet though.
        //if (!localize_first)
          throw runtime_error( string("Trouble assigning orbital to a monomer. |A|^2 = ") +
                               to_string(Anorm) + string(", |B|^2 = ") + to_string(Bnorm) );
      }
    }

    subsets_A.emplace_back(move(set_A));
    subsets_B.emplace_back(move(set_B));

    ambiguous_subsets.emplace_back(move(ambiguous));
  }

  auto out_coeff = local_coeff->copy();

  const int dimerbasis = sgeom_->nbasis();
  const int nsubspaces = orbital_subspaces.size();
  for (int sub = 0; sub < nsubspaces; ++sub) {
    size_t imo = orbital_subspaces[sub].first;

    vector<set<int>> subsets{{subsets_A[sub], subsets_B[sub], ambiguous_subsets[sub]}};
    for (auto& subset : subsets) {
      if (subset.empty()) continue;
      auto subspace = make_shared<Matrix>(dimerbasis, subset.size());
      int pos = 0;
      for (const int& i : subset)
        copy_n(local_coeff->element_ptr(0, i), dimerbasis, subspace->element_ptr(0, pos++));

      auto subfock = make_shared<Matrix>(*subspace % *fock * *subspace);
      vector<double> eigs(subspace->mdim());
      subfock->diagonalize(eigs.data());
      subspace = make_shared<Matrix>(*subspace * *subfock);

      copy_n(subspace->data(), dimerbasis * subset.size(), out_coeff->element_ptr(0,imo));
      imo += subset.size();
    }
  }

  sref_ = make_shared<Reference>(*sref_, make_shared<Coeff>(move(*out_coeff)));
}

/// localize_first flag defined as in Dimer::localize
///  If localize_first is set, separate SVDs are done for each monomer. Otherwise, two are done (one for
///  closed and one for active) on the whole dimer.
void Dimer::set_active(const std::shared_ptr<const PTree> idata, const bool localize_first) {
  // TODO needs clean up
  auto Asp = idata->get_child_optional("active_A");
  auto Bsp = idata->get_child_optional("active_B");
  auto isp = idata->get_child_optional("dimer_active");
  set<int> Ai, Bi, it;
  if (Asp) for (auto& s : *Asp) { Ai.insert(lexical_cast<int>(s->data())-1); } // TODO I think this -1 is very confusing!
  if (Bsp) for (auto& s : *Bsp) { Bi.insert(lexical_cast<int>(s->data())-1); }
  if (isp) for (auto& s : *isp) { it.insert(lexical_cast<int>(s->data())-1); }

  set<int> Alist, Blist;

  if (it.empty() && Ai.empty() && Bi.empty())
    throw runtime_error("Active space of the dimer MUST be specified in some way.");
  if (!it.empty()) {
    Alist = it;
    Blist = it;
  }
  if (!Ai.empty()) Alist = Ai;
  if (!Bi.empty()) Blist = Bi;

  // Make new References
  active_refs_ = make_pair(isolated_refs_.first->set_active(Alist), isolated_refs_.second->set_active(Blist));

  // Hold onto old occupation data
  const int noccA = isolated_refs_.first->nclosed();
  const int noccB = isolated_refs_.second->nclosed();

  const int nexternA = isolated_refs_.first->nvirt();
  const int nexternB = isolated_refs_.second->nvirt();

  // Update Dimer info
  const int nclosedA = active_refs_.first->nclosed();
  const int nclosedB = active_refs_.second->nclosed();
  const int nclosed = nclosedA + nclosedB;

  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;

  const int nactvirtA = isolated_refs_.first->nvirt() - active_refs_.first->nvirt();
  const int nactvirtB = isolated_refs_.second->nvirt() - active_refs_.second->nvirt();

  const int nbasisA = geoms_.first->nbasis();
  const int nbasisB = geoms_.second->nbasis();
  const int dimerbasis = sgeom_->nbasis();

  // TODO: this implementation requires specifying the number of active orbitals that are coming from each subset.
  //  This is probably fine, but it is not strictly necessary.

  // tuple:
  //  matrix --> reference active orbitals
  //  pair   --> bounds on subspace that can be mixed to produce localized orbitals
  //  int    --> number of orbitals to form
  //  string --> name of subspace (just for pretty printing)
  //  bool   --> closed(true)/virtual(false)
  vector<tuple<shared_ptr<const Matrix>, pair<int, int>, int, string, bool>> svd_info;

  if (localize_first) {
    auto activeA = make_shared<Matrix>(dimerbasis, nactA);
    activeA->copy_block(0, 0, nbasisA, nactA, active_refs_.first->coeff()->get_submatrix(0, nclosedA, nbasisA, nactA));
    svd_info.emplace_back(activeA, make_pair(0, noccA), noccA - nclosedA, "A", true);
    svd_info.emplace_back(activeA, make_pair(noccA+noccB, noccA+noccB+nexternA), nactvirtA, "A", false);

    auto activeB = make_shared<Matrix>(dimerbasis, nactB);
    activeB->copy_block(nbasisA, 0, nbasisB, nactB, active_refs_.second->coeff()->get_submatrix(0, nclosedB, nbasisB, nactB));
    svd_info.emplace_back(activeB, make_pair(noccA, noccA+noccB), noccB - nclosedB, "B", true);
    svd_info.emplace_back(activeB, make_pair(noccA+noccB+nexternA, noccA+noccB+nexternA+nexternB), nactvirtB, "B", false);
  }
  else {
    auto active = make_shared<Matrix>(dimerbasis, nact);

    active->copy_block(0, 0, nbasisA, nactA, active_refs_.first->coeff()->get_submatrix(0, nclosedA, nbasisA, nactA));
    active->copy_block(nbasisA, nactA, nbasisB, nactB, active_refs_.second->coeff()->get_submatrix(0, nclosedB, nbasisB, nactB));

    svd_info.emplace_back(active, make_pair(0, noccA + noccB), noccA + noccB - (nclosedA + nclosedB), "dimer", true);
    svd_info.emplace_back(active, make_pair(noccA + noccB, nbasisA + nbasisB), nactvirtA + nactvirtB, "dimer", false);
  }

  Overlap S(sgeom_);

  shared_ptr<Matrix> out_coeff = sref_->coeff()->copy();
  size_t closed_position = 0;
  size_t active_position = nclosed;
  size_t virt_position = nclosed + nact;

  for (auto& subset : svd_info) {
    const Matrix& active = *get<0>(subset);
    pair<int, int> bounds = get<1>(subset);
    const int norb = get<2>(subset);
    const string set_name = get<3>(subset);
    const bool closed = get<4>(subset);

    shared_ptr<Matrix> subcoeff = sref_->coeff()->slice_copy(bounds.first, bounds.second);

    const Matrix overlaps( active % S * *subcoeff );

    multimap<double, int> norms;

    for(int i = 0; i < overlaps.mdim(); ++i) {
      const double norm = blas::dot_product(overlaps.element_ptr(0, i), overlaps.ndim(), overlaps.element_ptr(0, i));
      norms.emplace(norm, i);
    }

    active_thresh_ = input_->get<double>("active_thresh", 0.5);
    cout << endl << "  o Forming dimer's active orbitals arising from " << (closed ? "closed " : "virtual ") <<  set_name << " orbitals. Threshold for inclusion in cadidate space: " << setw(6) << setprecision(3) << active_thresh_ << endl;

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
      Matrix subspace(dimerbasis, active_size);

      int ii = 0;
      for (int& i : active_list)
        copy_n(subcoeff->element_ptr(0, i), dimerbasis, subspace.element_ptr(0, ii++));

      Matrix Sactive(active % S * active);
      Sactive.inverse_half();

      Matrix projector( Sactive * ( active % S * subspace ) );
      vector<double> singulars(active_size, 0.0);
      shared_ptr<Matrix> Vt;
      tie(ignore, Vt) = projector.svd(singulars.data());

      cout << "    - largest singular value: " << singulars[0] << ", smallest: " << singulars[norb-1] << endl;
      cout << "    - norb: " << norb << ", sum of highest singular values: " << accumulate(singulars.begin(), singulars.begin()+norb, 0.0) << endl;

      subspace = subspace ^ *Vt;

      for (size_t i = 0; i < subcoeff->mdim(); ++i) {
        if ( count(active_list.begin(), active_list.end(), i) == 0 )
          copy_n(subcoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, ( closed ? closed_position++ : virt_position++ )));
      }
      copy_n(subspace.data(), dimerbasis * norb, out_coeff->element_ptr(0, active_position));
      active_position += norb;

      for (size_t i = norb; i < active_size; ++i)
        copy_n(subspace.element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, ( closed ? closed_position++ : virt_position++ )));
    }
    else {
      set<int> active_set(active_list.begin(), active_list.end());
      for (size_t i = 0; i < subcoeff->mdim(); ++i)
        if (active_set.count(i) == 0)
          copy_n(subcoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, ( closed ? closed_position++ : virt_position++ )));
        else
          copy_n(subcoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, active_position++));
    }
  }

  auto out = make_shared<Reference>(sgeom_, make_shared<Coeff>(*out_coeff), nclosed, nact, nexternA+nexternB - (nclosed+nact));


  sref_ = out;
}

// RHF and then localize
void Dimer::scf(const shared_ptr<const PTree> idata) {
  Timer dimertime;

  // SCF
  auto hfdata = idata->get_child_optional("hf") ? idata->get_child_optional("hf") : make_shared<PTree>();
  shared_ptr<SCF> rhf = dynamic_pointer_cast<SCF>(construct_method("hf", hfdata, sgeom_, sref_));
  rhf->compute();
  sref_ = rhf->conv_to_ref();
  dimertime.tick_print("Dimer SCF");

  const int nclosed = sref_->nclosed();

  shared_ptr<const Matrix> dimerdensity = sref_->coeff()->form_density_rhf(nclosed);
  shared_ptr<const Matrix> dimercoeff = sref_->coeff()->slice_copy(0,nclosed);

  // Explanation of schemes:
  //   localize_first           - fragment localizes, then picks the active space within each fragment (recommended)
  //   active_first             - picks active space from dimer orbitals first, then attempts to localize
  //   active_only              - picks active space from dimer orbitals and does not localize
  const string scheme = idata->get<string>("scheme", "active_first");

  if (scheme == "active_only" || scheme == "active_first") {
    // Set active space based on overlap
    if (isolated_refs_.first && isolated_refs_.second)
      set_active(idata, /*localize_first*/ false);
    else
      throw runtime_error("For Dimer::driver, Dimer must be constructed from HF references");

    if (scheme == "active_first") {
      shared_ptr<const PTree> localize_data = idata->get_child_optional("localization");
      if (!localize_data) localize_data = make_shared<const PTree>();

      auto fock  = make_shared<const Fock<1>>(sgeom_, sref_->hcore(), dimerdensity, dimercoeff);
      dimertime.tick_print("Dimer Fock matrix formation");

      localize(localize_data, fock, /*localize_first*/ false);
      dimertime.tick_print("Dimer localization");
    }
  }
  else if (scheme == "localize_first") {
    shared_ptr<const PTree> localize_data = idata->get_child_optional("localization");
    if (!localize_data) localize_data = make_shared<const PTree>();

    auto fock  = make_shared<const Fock<1>>(sgeom_, sref_->hcore(), dimerdensity, dimercoeff);
    dimertime.tick_print("Dimer Fock matrix formation");

    localize(localize_data, fock, /*localize_first*/ true);
    dimertime.tick_print("Dimer localization");

    set_active(idata, /*localize_first*/ true);

    const int nactA = active_refs_.first->nact();
    const int nactB = active_refs_.second->nact();
    Matrix active_mos = *sref_->coeff()->slice(nclosed, nclosed + nactA + nactB);
    Matrix fock_mo(active_mos % *fock * active_mos);
    vector<double> eigs(active_mos.mdim(), 0.0);
    shared_ptr<Matrix> active_transformation = fock_mo.diagonalize_blocks(eigs.data(), vector<int>{{nactA, nactB}});
    active_mos *= *active_transformation;
    shared_ptr<Matrix> scoeff = sref_->coeff()->copy();
    scoeff->copy_block(0, nclosed, scoeff->ndim(), active_mos.mdim(), active_mos);
    sref_ = make_shared<Reference>(*sref_, make_shared<Coeff>(move(*scoeff)));
  }
}


shared_ptr<DimerCAS> Dimer::compute_cispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  auto d1 = make_shared<Determinants>(active_refs_.first->nact(), nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = make_shared<Determinants>(active_refs_.second->nact(), nelea.second, neleb.second, /*compress*/false, /*mute*/true);
  auto out = make_shared<DimerCAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A;
  vector<vector<int>> spaces_B;

  auto space = idata->get_child_optional("space");
  if (space) {
    // TODO make a function
    for (auto& s : *space) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    // TODO make a function
    for (auto& s : *spacea) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    for (auto& s : *spaceb) { spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
  }

  Timer castime;

  shared_ptr<const PTree> fcidata = idata->get_child_optional("fci");
  if (!fcidata) fcidata = make_shared<const PTree>();

  // Embedded CAS-CI calculations
  cout << "    Starting embedded CAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_casci<0>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded CAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Charge, spin and number of states needs to be specified for each space");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_casci<1>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}


shared_ptr<DimerDistCAS> Dimer::compute_distcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  auto d1 = make_shared<Determinants>(active_refs_.first->nact(), nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = make_shared<Determinants>(active_refs_.second->nact(), nelea.second, neleb.second, /*compress*/false, /*mute*/true);
  auto out = make_shared<DimerDistCAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A;
  vector<vector<int>> spaces_B;

  auto space = idata->get_child_optional("space");
  if (space) {
    // TODO make a function
    for (auto& s : *space) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    // TODO make a function
    for (auto& s : *spacea) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    for (auto& s : *spaceb) { spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
  }

  Timer castime;

  shared_ptr<const PTree> fcidata = idata->get_child_optional("fci");
  if (!fcidata) fcidata = make_shared<const PTree>();

  // Embedded CAS-CI calculations
  cout << "    Starting embedded distributed CAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_distcasci<0>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded CAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Charge, spin and number of states needs to be specified for each space");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_distcasci<1>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}


shared_ptr<DimerRAS> Dimer::compute_rcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  // { {nras1, nras2, nras3}, max holes, max particles }
  pair<tuple<array<int, 3>, int, int>, tuple<array<int, 3>, int, int>> ras_desc;

  // Sample:
  // "restricted" : [ { "orbitals" : [1, 2, 3], "max_holes" : 0, "max_particles" : 2 } ],
  //  puts 1 orbital in RAS1 with no holes allowed, 2 orbital in RAS2, and 3 orbitals in RAS3 with up to 2 particles
  auto restrictions = idata->get_child("restricted");

  auto get_restricted_data = [] (shared_ptr<const PTree> i) {
    return make_tuple(i->get_array<int, 3>("orbitals"), i->get<int>("max_holes"), i->get<int>("max_particles"));
  };

  if (restrictions->size() == 1) {
    ras_desc = make_pair( get_restricted_data(*restrictions->begin()), get_restricted_data(*restrictions->begin()) );
  }
  else if (restrictions->size() == 2) {
    auto iter = restrictions->begin();
    auto tmp1 = get_restricted_data(*iter++);
    auto tmp2 = get_restricted_data(*iter);
    ras_desc = make_pair(tmp1, tmp2);
  }
  else throw logic_error("One or two sets of restrictions must be provided.");

  // This is less than ideal. It'd be better to have some sort of generator object that can be passed around.
  auto d1 = make_shared<RASDeterminants>(get<0>(ras_desc.first), nelea.first, neleb.first, get<1>(ras_desc.first), get<2>(ras_desc.first), true);
  auto d2 = make_shared<RASDeterminants>(get<0>(ras_desc.second), nelea.second, neleb.second, get<1>(ras_desc.second), get<2>(ras_desc.second), true);

  auto out = make_shared<DimerRAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A;
  vector<vector<int>> spaces_B;

  auto space = idata->get_child_optional("space");
  if (space) {
    // TODO make a function
    for (auto& s : *space) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    // TODO make a function
    for (auto& s : *spacea) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    for (auto& s : *spaceb) { spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
  }

  Timer castime;

  shared_ptr<const PTree> rasdata = idata->get_child_optional("ras");
  if (!rasdata) rasdata = make_shared<const PTree>();

  // Embedded RAS-CI calculations
  cout << "    Starting embedded RAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_rasci<0>(rasdata, charge, spin, nstate, ras_desc.first));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded RAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_rasci<1>(rasdata, charge, spin, nstate, ras_desc.second));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}


shared_ptr<DimerDistRAS> Dimer::compute_distrcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  // { {nras1, nras2, nras3}, max holes, max particles }
  pair<tuple<array<int, 3>, int, int>, tuple<array<int, 3>, int, int>> ras_desc;

  // Sample:
  // "restricted" : [ { "orbitals" : [1, 2, 3], "max_holes" : 0, "max_particles" : 2 } ],
  //  puts 1 orbital in RAS1 with no holes allowed, 2 orbital in RAS2, and 3 orbitals in RAS3 with up to 2 particles
  auto restrictions = idata->get_child("restricted");

  auto get_restricted_data = [] (shared_ptr<const PTree> i) {
    return make_tuple(i->get_array<int, 3>("orbitals"), i->get<int>("max_holes"), i->get<int>("max_particles"));
  };

  if (restrictions->size() == 1) {
    ras_desc = make_pair( get_restricted_data(*restrictions->begin()), get_restricted_data(*restrictions->begin()) );
  }
  else if (restrictions->size() == 2) {
    auto iter = restrictions->begin();
    auto tmp1 = get_restricted_data(*iter++);
    auto tmp2 = get_restricted_data(*iter);
    ras_desc = make_pair(tmp1, tmp2);
  }
  else throw logic_error("One or two sets of restrictions must be provided.");

  // This is less than ideal. It'd be better to have some sort of generator object that can be passed around.
  auto d1 = make_shared<RASDeterminants>(get<0>(ras_desc.first), nelea.first, neleb.first, get<1>(ras_desc.first), get<2>(ras_desc.first), true);
  auto d2 = make_shared<RASDeterminants>(get<0>(ras_desc.second), nelea.second, neleb.second, get<1>(ras_desc.second), get<2>(ras_desc.second), true);

  auto out = make_shared<DimerDistRAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A;
  vector<vector<int>> spaces_B;

  auto space = idata->get_child_optional("space");
  if (space) {
    // TODO make a function
    for (auto& s : *space) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    // TODO make a function
    for (auto& s : *spacea) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    for (auto& s : *spaceb) { spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
  }

  Timer castime;

  shared_ptr<const PTree> rasdata = idata->get_child_optional("ras");
  if (!rasdata) rasdata = make_shared<const PTree>();

  // Embedded RAS-CI calculations
  cout << "    Starting embedded RAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_distrasci<0>(rasdata, charge, spin, nstate, ras_desc.first));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded RAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_distrasci<1>(rasdata, charge, spin, nstate, ras_desc.second));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}
