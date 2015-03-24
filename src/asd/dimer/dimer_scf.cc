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

#include <src/asd/dimer/dimer.h>
#include <src/scf/hf/fock.h>
#include <src/wfn/localization.h>
#include <src/util/io/moldenout.h>

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

  pair<int, int> nbasis {geoms_.first->nbasis(), geoms_.second->nbasis()};

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
      VectorB eigs(subspace->mdim());
      subfock->diagonalize(eigs);
      subspace = make_shared<Matrix>(*subspace * *subfock);

      copy_n(subspace->data(), dimerbasis * subset.size(), out_coeff->element_ptr(0,imo));
      imo += subset.size();
    }
  }

  if (localize_first) {
    assert(nsubspaces==2);
    if (orbital_subspaces[0].first > orbital_subspaces[1].first)
      nvirt_ = {subsets_A[0].size(), subsets_B[0].size()};
    else
      nvirt_ = {subsets_A[1].size(), subsets_B[1].size()};
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
  active_refs_ = {isolated_refs_.first->set_active(Alist), isolated_refs_.second->set_active(Blist)};

  // Hold onto old occupation data
  const int noccA = isolated_refs_.first->nclosed();
  const int noccB = isolated_refs_.second->nclosed();

  const int nexternA = nvirt_.first;
  const int nexternB = nvirt_.second;

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

  nvirt_ = {nexternA - nactA, nexternB - nactB};
  sref_ = make_shared<Reference>(sgeom_, make_shared<Coeff>(*out_coeff), nclosed, nact, nexternA+nexternB - (nclosed+nact));
}


// RHF and then localize
void Dimer::scf(const shared_ptr<const PTree> idata) {
  Timer dimertime;

  // SCF
  auto hfdata = idata->get_child_optional("hf") ? idata->get_child_optional("hf") : make_shared<PTree>();
  auto rhf = dynamic_pointer_cast<RHF>(construct_method("hf", hfdata, sgeom_, sref_));
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
  //   linked                   - similar to localize_first, but for linked dimer only
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
    Matrix active_mos = sref_->coeff()->slice(nclosed, nclosed + nactA + nactB);
    Matrix fock_mo(active_mos % *fock * active_mos);
    VectorB eigs(active_mos.mdim());
    shared_ptr<Matrix> active_transformation = fock_mo.diagonalize_blocks(eigs, vector<int>{{nactA, nactB}});
    active_mos *= *active_transformation;
    shared_ptr<Matrix> scoeff = sref_->coeff()->copy();
    scoeff->copy_block(0, nclosed, scoeff->ndim(), active_mos.mdim(), active_mos);
    sref_ = make_shared<Reference>(*sref_, make_shared<Coeff>(move(*scoeff)));
  }
  else if (scheme == "legacy_mode") { //TODO: to be removed, and replaced by "linked"
    assert(idata->get<string>("form") == "legacy_mode");
    shared_ptr<const PTree> localize_data = idata->get_child_optional("localization");
    if (!localize_data) localize_data = make_shared<const PTree>();

    auto fock  = make_shared<const Fock<1>>(sgeom_, sref_->hcore(), dimerdensity, dimercoeff);
    dimertime.tick_print("Dimer Fock matrix formation");

    string localizemethod = localize_data->get<string>("algorithm", "pm");
    shared_ptr<OrbitalLocalization> localization;
    if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
      localization = make_shared<PMLocalization>(localize_data, sref_);
    else throw runtime_error("Unrecognized orbital localization method");

    shared_ptr<const Coeff> new_coeff = make_shared<const Coeff>(*localization->localize());
    sref_ = make_shared<const Reference>(*sref_, new_coeff);

    shared_ptr<const PTree> printAB_data = idata->get_child_optional("print_AB");
    if (printAB_data) {
      const bool orbitals = printAB_data->get<bool>("orbitals", true);
      const string out_file = printAB_data->get<string>("file", "AB.molden");

      if (mpi__->rank() == 0) {
        MoldenOut mfs(out_file);
        mfs << sgeom_;
        if (orbitals) mfs << sref_;
      }
    }

    auto Asp = idata->get_child_optional("active_A");
    auto Bsp = idata->get_child_optional("active_B");
    set<int> Alist, Blist;
    if (Asp) for (auto& s : *Asp) { Alist.insert(lexical_cast<int>(s->data())-1); }
    if (Bsp) for (auto& s : *Bsp) { Blist.insert(lexical_cast<int>(s->data())-1); }

    if (Alist.empty() && Blist.empty())
      throw runtime_error("Active space of the dimer MUST be specified in some way.");

    const int nactA = Alist.size();
    const int nactB = Blist.size();
    const int nbasis = sgeom_->nbasis();

    set<int> ABlist;
    if (!Alist.empty()) for (auto i : Alist) { ABlist.insert(i); }
    if (!Blist.empty()) for (auto i : Blist) { ABlist.insert(i); }
    assert(nactA+nactB == ABlist.size());

    isolated_refs_ = {make_shared<Reference>(*sref_), make_shared<Reference>(*sref_)}; //sref_: localized HF MO's, not sorted

    active_refs_ = {isolated_refs_.first->set_active(Alist), isolated_refs_.second->set_active(Blist)};

    shared_ptr<const PTree> printA_data = idata->get_child_optional("print_A");
    if (printA_data) {
      const bool orbitals = printA_data->get<bool>("orbitals", true);
      const string out_file = printA_data->get<string>("file", "A.molden");

      if (mpi__->rank() == 0) {
        MoldenOut mfs(out_file);
        mfs << sgeom_;
        if (orbitals) mfs << active_refs_.first;
      }
    }

    shared_ptr<const PTree> printB_data = idata->get_child_optional("print_B");
    if (printB_data) {
      const bool orbitals = printB_data->get<bool>("orbitals", true);
      const string out_file = printB_data->get<string>("file", "B.molden");

      if (mpi__->rank() == 0) {
        MoldenOut mfs(out_file);
        mfs << sgeom_;
        if (orbitals) mfs << active_refs_.second;
      }
    }

    //activate sref_: order = (nclosed, nactA, nactB, nvirt)
    auto tmpref = sref_->set_active(ABlist);
    auto tmpcoeff = make_shared<Matrix>(*tmpref->coeff());
    shared_ptr<const Matrix> coeffA = active_refs_.first->coeff();
    shared_ptr<const Matrix> coeffB = active_refs_.second->coeff();
    const int nclosed = tmpref->nclosed();
    const int nclosedA = active_refs_.first->nclosed();
    const int nclosedB = active_refs_.second->nclosed();
    tmpcoeff->copy_block(0, nclosed, nbasis, nactA, coeffA->slice(nclosedA,nclosedA+nactA));
    tmpcoeff->copy_block(0, nclosed+nactA, nbasis, nactB, coeffB->slice(nclosedB,nclosedB+nactB));
    sref_ = make_shared<Reference>(*tmpref, make_shared<Coeff>(move(*tmpcoeff)));
  }
  else if (scheme == "linked") {
    shared_ptr<const PTree> localize_data = idata->get_child_optional("localization");
    if (!localize_data) localize_data = make_shared<const PTree>();

    string localizemethod = localize_data->get<string>("algorithm", "pm");
    shared_ptr<OrbitalLocalization> localization;
    if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
      localization = make_shared<PMLocalization>(localize_data, sref_);
    else throw runtime_error("Unrecognized orbital localization method");

    shared_ptr<const Coeff> new_coeff = make_shared<const Coeff>(*localization->localize());
    sref_ = make_shared<const Reference>(*sref_, new_coeff); //new super_reference with localized orbitals

#if 1
    if (mpi__->rank() == 0) {
      MoldenOut mfs("AB_localized_L.molden");
      mfs << sgeom_;
      mfs << sref_;
    }
#endif

    set_active(idata, /*not used here*/true, /*linked*/true);

  }
}


//TODO: merge
void Dimer::set_active(const std::shared_ptr<const PTree> idata, const bool localize_first, const bool linked) {
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
    vector<pair<int, int>> bounds;
    set<int> ambiguous_list;
    //sort according to region and put into Alist & Blist
    vector<int> sizes = idata->get_vector<int>("region_sizes");
    int nbasis = 0;
    int natoms = 0;
    for (int& region : sizes) {
      const int atomstart = natoms;
      const int basisstart = nbasis;
      for (int atom = atomstart; atom < atomstart + region; ++atom)
        nbasis += sgeom_->atoms()[atom]->nbasis();

      natoms += region;
      if (basisstart != nbasis)
        bounds.emplace_back(basisstart, nbasis);
    }
    if (bounds.size() != 3)
      throw logic_error("Only 3 regions should be defined in order of A, B, and the rest");
    if (natoms != count_if(sgeom_->atoms().begin(), sgeom_->atoms().end(), [](const shared_ptr<const Atom> a){return !a->dummy();}))
      throw logic_error("All atoms must be assigned to regions");
    //scan it set
    cout << "  o Assigning dimer active orbitals (localized) to monomers A and B" << endl;
    shared_ptr<Matrix> coeff = isolated_refs_.first->coeff()->copy(); //projected coeff on isolated_refs (so far, first == second respresenting dimer hf)
    for (auto& imo : it) {
      const double sum_A = blas::dot_product(coeff->element_ptr(bounds[0].first, imo), bounds[0].second - bounds[0].first, coeff->element_ptr(bounds[0].first, imo));
      const double sum_B = blas::dot_product(coeff->element_ptr(bounds[1].first, imo), bounds[1].second - bounds[1].first, coeff->element_ptr(bounds[1].first, imo));
      const double sum_rest = blas::dot_product(coeff->element_ptr(bounds[2].first, imo), bounds[2].second - bounds[2].first, coeff->element_ptr(bounds[2].first, imo));
      if (sum_A > sum_B && sum_A > sum_rest) {
        cout << "    - orbital("  << imo << ") is assigned to monomer A : A(" << setw(6) << setprecision(3) << sum_A << "), B("
                                                                              << setw(6) << setprecision(3) << sum_B << "), the rest("
                                                                              << setw(6) << setprecision(3) << sum_rest << ")" << endl;
        Alist.insert(imo);
      }
      else if (sum_B > sum_A && sum_B > sum_rest) {
        cout << "    - orbital("  << imo << ") is assigned to monomer B : A(" << setw(6) << setprecision(3) << sum_A << "), B("
                                                                              << setw(6) << setprecision(3) << sum_B << "), the rest("
                                                                              << setw(6) << setprecision(3) << sum_rest << ")" << endl;
        Blist.insert(imo);
      }
      else {
        cout << "    - orbital("  << imo << ") cannot be assigned to neither of monomers : A(" << setw(6) << setprecision(3) << sum_A << "), B("
                                                                              << setw(6) << setprecision(3) << sum_B << "), the rest("
                                                                              << setw(6) << setprecision(3) << sum_rest << ")" << endl;
        throw runtime_error("Wrong choice of active orbitals");
      }
    }
  }
  if (!Ai.empty()) Alist = Ai;
  if (!Bi.empty()) Blist = Bi;

  // Make new References, with large basis sets, but with projected coeffs, active orbitals are placed after closed orbitals
  active_refs_ = {isolated_refs_.first->set_active(Alist), isolated_refs_.second->set_active(Blist)};

#if 1
  if (mpi__->rank() == 0) {
    MoldenOut mfs("A_active_L.molden");
    mfs << active_refs_.first->geom();
    mfs << active_refs_.first;
  }
  if (mpi__->rank() == 0) {
    MoldenOut mfs("B_active_L.molden");
    mfs << active_refs_.second->geom();
    mfs << active_refs_.second;
  }
#endif

  // Update Dimer info
  const int nclosedA = active_refs_.first->nclosed();
  const int nclosedB = active_refs_.second->nclosed(); //these shares common closed indices, but with different embedded active(closed part) orbtals

  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;

  const int nactvirtA = isolated_refs_.first->nvirt() - active_refs_.first->nvirt();
  const int nactvirtB = isolated_refs_.second->nvirt() - active_refs_.second->nvirt(); // Number of active orbitals in virtual HF subspace

  const int dimerbasis = sgeom_->nbasis();
  assert(dimerbasis == geoms_.first->nbasis());
  assert(dimerbasis == geoms_.second->nbasis());

  const int nclosedS = sref_->nclosed();
  const int nvirtS = sref_->nvirt();
  assert(dimerbasis == nclosedS + nvirtS);
  assert(sref_->nact() == 0);

  const int nclosed = nclosedS - (nclosedS - nclosedA) - (nclosedS - nclosedB);

  // TODO: this implementation requires specifying the number of active orbitals that are coming from each subset.
  //  This is probably fine, but it is not strictly necessary.

  // tuple:
  //  matrix --> reference active orbitals
  //  pair   --> bounds on subspace that can be mixed to produce localized orbitals
  //  int    --> number of orbitals to form
  //  string --> name of subspace (just for pretty printing)
  //  bool   --> closed(true)/virtual(false)
  vector<tuple<shared_ptr<const Matrix>, pair<int, int>, int, string, bool>> svd_info;

  auto activeA = make_shared<Matrix>(dimerbasis, nactA);
  activeA->copy_block(0, 0, dimerbasis, nactA, active_refs_.first->coeff()->get_submatrix(0, nclosedA, dimerbasis, nactA));
  svd_info.emplace_back(activeA, make_pair(0, nclosedS), nclosedS - nclosedA, "A", true); //A active in closed subspace
  svd_info.emplace_back(activeA, make_pair(nclosedS, dimerbasis), nactvirtA, "A", false); //A active in virtual subspace

  auto activeB = make_shared<Matrix>(dimerbasis, nactB);
  activeB->copy_block(0, 0, dimerbasis, nactB, active_refs_.second->coeff()->get_submatrix(0, nclosedB, dimerbasis, nactB));
  svd_info.emplace_back(activeB, make_pair(0, nclosedS), nclosedS - nclosedB, "B", true);
  svd_info.emplace_back(activeB, make_pair(nclosedS, dimerbasis), nactvirtB, "B", false);

  Overlap S(sgeom_);

  shared_ptr<Matrix> out_coeff = sref_->coeff()->copy();
  size_t active_position = nclosed;

  set<int> mask; //records what orbitals have been copied to sref_ coeff, to be used to copy back the common closed & virtual orbitals
  for (int i = 0; i != out_coeff->mdim(); ++i) mask.insert(i);

  for (auto& subset : svd_info) {
    const Matrix& active = *get<0>(subset);
    pair<int, int> bounds = get<1>(subset);
    const int norb = get<2>(subset);
    const string set_name = get<3>(subset);
    const bool closed = get<4>(subset);

    shared_ptr<Matrix> subcoeff = sref_->coeff()->slice_copy(bounds.first, bounds.second);

    const Matrix overlaps(active % S * *subcoeff);

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
      throw runtime_error("SVD is not yet implemented. Try adjust active_thresh.");
    }
    else {
      set<int> active_set(active_list.begin(), active_list.end());
      for (size_t i = 0; i < subcoeff->mdim(); ++i) {
        if (active_set.count(i)) {
          const int imo = bounds.first + i;
          assert(mask.count(imo));
          mask.erase(imo);
          copy_n(subcoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, active_position++));
        }
      }
    }
  }

  //fill common closed and virtual subspace
  shared_ptr<Matrix> scoeff = sref_->coeff()->copy();

  size_t closed_position = 0;
  for (int i = 0; i < nclosedS; ++i)
    if(mask.count(i))
      copy_n(scoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, closed_position++));

  size_t virt_position = nclosed + nact;
  for (int i = nclosedS; i < dimerbasis; ++i)
    if(mask.count(i))
      copy_n(scoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, virt_position++));

  nvirt_ = {nvirtS - nactvirtA, nvirtS - nactvirtB};
  sref_ = make_shared<Reference>(sgeom_, make_shared<Coeff>(*out_coeff), nclosed, nact, nvirtS - nactvirtA - nactvirtB);

  if (mpi__->rank() == 0) {
    MoldenOut mfs("AB_localized_reordered_L.molden");
    mfs << sgeom_;
    mfs << sref_;
  }
}
