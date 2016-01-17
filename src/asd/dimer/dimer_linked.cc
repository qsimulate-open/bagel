//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/dimer/dimer_linked.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/dimer/dimer.h>
#include <src/scf/hf/fock.h>
#include <src/wfn/localization.h>
#include <src/util/io/moldenout.h>

using namespace std;
using namespace bagel;

void Dimer::set_active(shared_ptr<const PTree> idata) {
  auto isp = idata->get_child_optional("dimer_active");
  set<int> it;
  if (isp) for (auto& s : *isp) { it.insert(lexical_cast<int>(s->data())-1); }
  if (it.empty())
    throw runtime_error("Active space of the dimer MUST be specified in dimer_active.");

  set<int> Alist, Blist;

  { //sort out mixed active indices into A or B fragment (defined by region), if an index does not belong to any of fragments, throw error..
    vector<pair<int, int>> bounds;
    vector<int> sizes = idata->get_vector<int>("region_sizes"); // [A, B, the rest]
    int nbasis = 0;
    int natoms = 0;
    //set bounds : first & last basis function indices for the given region
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
    cout << "  o Assigning dimer active (localized & projected) orbitals to monomers A and B" << endl;
    shared_ptr<Matrix> coeff = isolated_refs_.first->coeff()->copy(); //"projected coeff" on isolated_refs (so far, first == second respresenting dimer hf)
    for (auto& imo : it) {
      const double sum_A = blas::dot_product(coeff->element_ptr(bounds[0].first, imo), bounds[0].second - bounds[0].first, coeff->element_ptr(bounds[0].first, imo));
      const double sum_B = blas::dot_product(coeff->element_ptr(bounds[1].first, imo), bounds[1].second - bounds[1].first, coeff->element_ptr(bounds[1].first, imo));
      const double sum_rest = blas::dot_product(coeff->element_ptr(bounds[2].first, imo), bounds[2].second - bounds[2].first, coeff->element_ptr(bounds[2].first, imo));

      if (sum_A > sum_B && sum_A > sum_rest)
      //cout << "    - active orbital("  << imo << ") is assigned to monomer A." << endl;
      //cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << "), the rest(" << setw(6) << setprecision(3) << sum_rest << ")" << endl;
        Alist.insert(imo);
      else if (sum_B > sum_A && sum_B > sum_rest)
        Blist.insert(imo);
      else
        throw runtime_error("Wrong choice of active orbitals. The orbital(" + to_string(imo) + ") does not belong to any of monomers.");
    }
    cout << "    - orbitals are assigned as : " << Alist.size() << "(A), " << Blist.size() << "(B)." << endl;
  }

  // Make new References, with large basis sets, but with projected coeffs (MO space up to smaller basis sets); active orbitals are placed after closed orbitals
  active_refs_ = {isolated_refs_.first->set_active(Alist), isolated_refs_.second->set_active(Blist)};

  if (print_orbital_ && mpi__->rank() == 0) {
    {
      MoldenOut mfs("dimer_orbital_A_projected.molden");
      mfs << active_refs_.first->geom();
      mfs << active_refs_.first;
    }
    {
      MoldenOut mfs("dimer_orbital_B_projected.molden");
      mfs << active_refs_.second->geom();
      mfs << active_refs_.second;
    }
  }

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
  //so far, sref_ is dimer RHF reference
  const int nclosed_HF = sref_->nclosed();
  const int nvirt_HF = sref_->nvirt();
  assert(dimerbasis == nclosed_HF + nvirt_HF);
  assert(sref_->nact() == 0);
  const int nclosed = nclosed_HF - (nclosed_HF - nclosedA) - (nclosed_HF - nclosedB);

  //prepare reference active orbitals
  auto control = make_shared<Matrix>(dimerbasis,nact);
  control->copy_block(0,0,     dimerbasis,nactA, active_refs_.first->coeff()->get_submatrix(0,nclosedA,dimerbasis,nactA));
  control->copy_block(0,nactA, dimerbasis,nactB, active_refs_.second->coeff()->get_submatrix(0,nclosedB,dimerbasis,nactB));

  //pick active orbitals based on reference and reorder to [closed,A,B,virt]
  cout << endl << "  o Picking up active orbitals using projected orbitals" << endl;
  shared_ptr<Matrix> out_coeff = overlap_selection(control, sref_->coeff());
  nvirt_ = {nvirt_HF - nactvirtA, nvirt_HF - nactvirtB};
  sref_ = make_shared<Reference>(sgeom_, make_shared<Coeff>(*out_coeff), nclosed, nact, nvirt_HF - nactvirtA - nactvirtB);
}


//semi-canonicalise only in the active space (thus, it is still the same quality as the localized orbitals)
//this will be used as reference to find the actual semi-canonical orbital
shared_ptr<Matrix> Dimer::form_reference_active_coeff() const {
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;
  const int nactcloA = isolated_refs_.first->nclosed()  - active_refs_.first->nclosed();
  const int nactcloB = isolated_refs_.second->nclosed() - active_refs_.second->nclosed();
  const int nclosed = sref_->nclosed();
  const int dimerbasis = sgeom_->nbasis();

  auto active_semi_coeff = make_shared<Matrix>(dimerbasis,nact);

  //closed
  auto ccoeff = make_shared<Matrix>(dimerbasis, nclosed+nactcloA+nactcloB);
  ccoeff->copy_block(0,0,                dimerbasis,nclosed,  sref_->coeff()->get_submatrix(0,0,             dimerbasis,nclosed)); //shared closed
  ccoeff->copy_block(0,nclosed,          dimerbasis,nactcloA, sref_->coeff()->get_submatrix(0,nclosed,       dimerbasis,nactcloA)); //embed activeB
  ccoeff->copy_block(0,nclosed+nactcloA, dimerbasis,nactcloB, sref_->coeff()->get_submatrix(0,nclosed+nactA, dimerbasis,nactcloB)); //embed activeB

  //AO Fock
  shared_ptr<const Matrix> ofockao = make_shared<Fock<1>>(sgeom_, sref_->hcore(), nullptr, ccoeff, /*store*/false, /*rhf*/true);

  {//Monomer A
    auto acoeff = sref_->coeff()->slice_copy(nclosed, nclosed+nactA);
    // MO Fock
    VectorB eigs(nactA);
    auto fockact = make_shared<Matrix>(*acoeff % *ofockao  * *acoeff);
    fockact->diagonalize(eigs);
    cout << endl << "  o Eigenvlues of A orbitals :" << endl;
    for (int i = 0; i < nactA; ++i) cout << setw(12) << setprecision(6) << eigs[i];
    cout << endl << endl;
    *acoeff *= *fockact;

    size_t act_position = 0; //for A
    for (int i = 0; i < nactA; ++i)
    copy_n(acoeff->element_ptr(0, i), dimerbasis, active_semi_coeff->element_ptr(0,act_position++));
  }

  {//Monomer B
    auto acoeff = sref_->coeff()->slice_copy(nclosed+nactA, nclosed+nact);
    // MO Fock
    VectorB eigs(nactB);
    auto fockact = make_shared<Matrix>(*acoeff % *ofockao  * *acoeff);
    fockact->diagonalize(eigs);
    cout << "  o Eigenvlues of B orbitals :" << endl;
    for (int i = 0; i < nactB; ++i) cout << setw(12) << setprecision(6) << eigs[i];
    cout << endl << endl;
    *acoeff *= *fockact;

    size_t act_position = nactA; //for B
    for (int i = 0; i < nactB; ++i)
      copy_n(acoeff->element_ptr(0, i), dimerbasis, active_semi_coeff->element_ptr(0,act_position++));
  }

  return active_semi_coeff;
}


shared_ptr<Matrix> Dimer::form_semi_canonical_coeff(shared_ptr<const PTree> idata) const {
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;
  const int nactcloA = isolated_refs_.first->nclosed()  - active_refs_.first->nclosed();
  const int nactcloB = isolated_refs_.second->nclosed() - active_refs_.second->nclosed();
  const int nactvirtA = isolated_refs_.first->nvirt() - active_refs_.first->nvirt();
  const int nactvirtB = isolated_refs_.second->nvirt() - active_refs_.second->nvirt(); // Number of active orbitals in virtual HF subspace
  const int nclosed = sref_->nclosed();
  const int nclosed_HF = sref_->nclosed() + nactcloA + nactcloB;
  const int nvirt_HF = sref_->nvirt();
  const int dimerbasis = sgeom_->nbasis();

  //completely semi-canonicalized based on fragments. [ closed | virtual]
  shared_ptr<Matrix> semi_coeff = sref_->coeff()->copy();

  //sort closed & virtual into A, B and bridge sites
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
  set<int> closed_Alist, closed_Blist, closed_Clist;
  set<int> virtual_Alist, virtual_Blist, virtual_Clist;
  shared_ptr<Matrix> coeff = semi_coeff->copy();
  {//closed
    set<int> cset; // (0:nclosed)
    for (int i = 0; i != nclosed; ++i) cset.insert(i);
    cout << "  o Assigning dimer closed (localized) orbitals to monomers A and B, and bridge" << endl;
    for (auto& imo : cset) {
      const double sum_A = blas::dot_product(coeff->element_ptr(bounds[0].first, imo), bounds[0].second - bounds[0].first, coeff->element_ptr(bounds[0].first, imo));
      const double sum_B = blas::dot_product(coeff->element_ptr(bounds[1].first, imo), bounds[1].second - bounds[1].first, coeff->element_ptr(bounds[1].first, imo));
      const double sum_rest = blas::dot_product(coeff->element_ptr(bounds[2].first, imo), bounds[2].second - bounds[2].first, coeff->element_ptr(bounds[2].first, imo));
      if (sum_A > sum_B && sum_A > sum_rest)
        closed_Alist.insert(imo);
      //cout << "    - closed orbital("  << imo << ") is assigned to monomer A." << endl;
      //cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << "), the rest(" << setw(6) << setprecision(3) << sum_rest << ")" << endl;
      else if (sum_B > sum_A && sum_B > sum_rest)
        closed_Blist.insert(imo);
      else
        closed_Clist.insert(imo);
    }
    cout << "    - orbitals are assigned as : " << closed_Alist.size() << "(A), " << closed_Blist.size() << "(B), " << closed_Clist.size() << "(bridge) " << endl << endl;
  }
  {//virtual
    set<int> vset; // (nclosed+nact:nbasis)
    for (int i = nclosed+nact; i != nbasis; ++i) vset.insert(i);
    cout << "  o Assigning dimer virtual (localized) orbitals to monomers A and B, and bridge" << endl;
    for (auto& imo : vset) {
      const double sum_A = blas::dot_product(coeff->element_ptr(bounds[0].first, imo), bounds[0].second - bounds[0].first, coeff->element_ptr(bounds[0].first, imo));
      const double sum_B = blas::dot_product(coeff->element_ptr(bounds[1].first, imo), bounds[1].second - bounds[1].first, coeff->element_ptr(bounds[1].first, imo));
      const double sum_rest = blas::dot_product(coeff->element_ptr(bounds[2].first, imo), bounds[2].second - bounds[2].first, coeff->element_ptr(bounds[2].first, imo));
      if (sum_A > sum_B && sum_A > sum_rest)
        virtual_Alist.insert(imo);
      //cout << "    - virtual orbital("  << imo << ") is assigned to monomer A." << endl;
      else if (sum_B > sum_A && sum_B > sum_rest)
        virtual_Blist.insert(imo);
      else
        virtual_Clist.insert(imo);
    }
    cout << "    - orbitals are assigned as : " << virtual_Alist.size() << "(A), " << virtual_Blist.size() << "(B), " << virtual_Clist.size() << "(bridge) " << endl;
  }

  auto ccoeff_A = make_shared<Matrix>(dimerbasis,nclosed_HF);
  auto ccoeff_B = ccoeff_A->clone();
  auto ccoeff_C = ccoeff_A->clone();

  auto vcoeff_A = make_shared<Matrix>(dimerbasis,nvirt_HF);
  auto vcoeff_B = vcoeff_A->clone();
  auto vcoeff_C = vcoeff_A->clone();

  //Put orbitals to appropriate coefficient matrices
  //closed + active
  {//Monomer A
    size_t pos = 0;
    for (auto& i : closed_Alist)
      copy_n(coeff->element_ptr(0, i), dimerbasis, ccoeff_A->element_ptr(0,pos++));
    for (int i = 0; i != nactcloA; ++i)
      copy_n(coeff->element_ptr(0, nclosed+i), dimerbasis, ccoeff_A->element_ptr(0,pos++));
  }
  {//Monomer B
    size_t pos = 0;
    for (auto& i : closed_Blist)
      copy_n(coeff->element_ptr(0, i), dimerbasis, ccoeff_B->element_ptr(0,pos++));
    for (int i = 0; i != nactcloB; ++i)
      copy_n(coeff->element_ptr(0, nclosed+nactA+i), dimerbasis, ccoeff_B->element_ptr(0,pos++));
  }
  {//Bridge
    size_t pos = 0;
    for (auto& i : closed_Clist)
      copy_n(coeff->element_ptr(0, i), dimerbasis, ccoeff_C->element_ptr(0,pos++));
  }
  //virtual + active
  {//Monomer A
    size_t pos = 0;
    for (auto& i : virtual_Alist)
      copy_n(coeff->element_ptr(0, i), dimerbasis, vcoeff_A->element_ptr(0,pos++));
    for (int i = 0; i != nactvirtA; ++i)
      copy_n(coeff->element_ptr(0, nclosed+nactcloA+i), dimerbasis, vcoeff_A->element_ptr(0,pos++));
  }
  {//Monomer B
    size_t pos = 0;
    for (auto& i : virtual_Blist)
      copy_n(coeff->element_ptr(0, i), dimerbasis, vcoeff_B->element_ptr(0,pos++));
    for (int i = 0; i != nactvirtB; ++i)
      copy_n(coeff->element_ptr(0, nclosed+nactA+nactcloB+i), dimerbasis, vcoeff_B->element_ptr(0,pos++));
  }
  {//Bridge
    size_t pos = 0;
    for (auto& i : virtual_Clist)
      copy_n(coeff->element_ptr(0, i), dimerbasis, vcoeff_C->element_ptr(0,pos++));
  }

  //form AO Fock
  shared_ptr<const Matrix> ofockao;
  {
    assert(nclosed_HF == nclosed + nactcloA + nactcloB);
    auto ccoeff = make_shared<Matrix>(dimerbasis, nclosed_HF);
    ccoeff->copy_block(0,0,                dimerbasis,nclosed,   *coeff->get_submatrix(0,0,             dimerbasis,nclosed)); //shared closed
    ccoeff->copy_block(0,nclosed,          dimerbasis,nactcloA,  *coeff->get_submatrix(0,nclosed,       dimerbasis,nactcloA)); //embed activeA
    ccoeff->copy_block(0,nclosed+nactcloA, dimerbasis,nactcloB,  *coeff->get_submatrix(0,nclosed+nactA, dimerbasis,nactcloB)); //embed activeB
    ofockao = make_shared<Fock<1>>(sgeom_, sref_->hcore(), nullptr, ccoeff, /*store*/false, /*rhf*/true);
  }

  //form MO fock
  size_t pos = 0; // order=[closed(A,B,C)|virtual(A,B,C)]
  { //closed A
    const int dim = closed_Alist.size()+nactcloA;
    VectorB eigs(dim);
    auto mocoeff = ccoeff_A->slice_copy(0,dim);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao  * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock; //transformed
    for (int i = 0; i < dim; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, semi_coeff->element_ptr(0,pos++)); //put back to right place of output matrix
  }
  { //closed B
    const int dim = closed_Blist.size()+nactcloB;
    VectorB eigs(dim);
    auto mocoeff = ccoeff_B->slice_copy(0,dim);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao  * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock; //transformed
    for (int i = 0; i < dim; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, semi_coeff->element_ptr(0,pos++));
  }
  if (!closed_Clist.empty()) { //closed C
    const int dim = closed_Clist.size();
    VectorB eigs(dim);
    auto mocoeff = ccoeff_C->slice_copy(0,dim);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao  * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock; //transformed
    for (int i = 0; i < dim; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, semi_coeff->element_ptr(0,pos++));
  }
  { //virtual A
    const int dim = virtual_Alist.size()+nactvirtA;
    VectorB eigs(dim);
    auto mocoeff = vcoeff_A->slice_copy(0,dim);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao  * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock; //transformed
    for (int i = 0; i < dim; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, semi_coeff->element_ptr(0,pos++));
  }
  { //virtual B
    const int dim = virtual_Blist.size()+nactvirtB;
    VectorB eigs(dim);
    auto mocoeff = vcoeff_B->slice_copy(0,dim);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao  * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock; //transformed
    for (int i = 0; i < dim; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, semi_coeff->element_ptr(0,pos++));
  }
  if (!virtual_Clist.empty()) { //virtual C
    const int dim = virtual_Clist.size();
    VectorB eigs(dim);
    auto mocoeff = vcoeff_C->slice_copy(0,dim);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao  * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock; //transformed
    for (int i = 0; i < dim; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, semi_coeff->element_ptr(0,pos++));
  }

  return semi_coeff;
}


shared_ptr<Matrix> Dimer::overlap_selection(shared_ptr<const Matrix> control, shared_ptr<const Matrix> treatment) const {
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;

  const int nactcloA = isolated_refs_.first->nclosed()  - active_refs_.first->nclosed();
  const int nactcloB = isolated_refs_.second->nclosed() - active_refs_.second->nclosed();
  const int nactvirtA = isolated_refs_.first->nvirt() - active_refs_.first->nvirt();
  const int nactvirtB = isolated_refs_.second->nvirt() - active_refs_.second->nvirt(); // Number of active orbitals in virtual HF subspace
  assert(nactA == nactcloA + nactvirtA);
  assert(nactB == nactcloB + nactvirtB);

  const int dimerbasis = sgeom_->nbasis();
  const int nclosed_HF = isolated_refs_.first->nclosed();
  const int nclosed = nclosed_HF - nactcloA - nactcloB;

  vector<tuple<shared_ptr<const Matrix>, pair<int, int>, int, string, bool>> ovl_info;

  auto activeA = make_shared<Matrix>(dimerbasis, nactA);
  activeA->copy_block(0, 0, dimerbasis, nactA, control->get_submatrix(0,0, dimerbasis, nactA));
  ovl_info.emplace_back(activeA, make_pair(0, nclosed_HF), nactcloA, "A", true); //A active in closed subspace
  ovl_info.emplace_back(activeA, make_pair(nclosed_HF, dimerbasis), nactvirtA, "A", false); //A active in virtual subspace

  auto activeB = make_shared<Matrix>(dimerbasis, nactB);
  activeB->copy_block(0, 0, dimerbasis, nactB, control->get_submatrix(0,nactA, dimerbasis, nactB));
  ovl_info.emplace_back(activeB, make_pair(0, nclosed_HF), nactcloB, "B", true);
  ovl_info.emplace_back(activeB, make_pair(nclosed_HF, dimerbasis), nactvirtB, "B", false);

  const Overlap S(sgeom_);

  shared_ptr<Matrix> out_coeff = treatment->clone();
  size_t active_position = nclosed;

  set<int> mask; //records what orbitals have been copied to sref_ coeff, to be used to copy back the common closed & virtual orbitals
  for (int i = 0; i != out_coeff->mdim(); ++i) mask.insert(i);

  //this fills the active (A and B) in order of closed A - virtual A - closed B - virtual B
  for (auto& subset : ovl_info) {
    const Matrix& active = *get<0>(subset);
    const pair<int, int> bounds = get<1>(subset);
    const int norb = get<2>(subset);
    const string set_name = get<3>(subset);
    const bool closed = get<4>(subset);

    shared_ptr<const Matrix> subcoeff = treatment->slice_copy(bounds.first, bounds.second);

    const Matrix overlaps(active % S * *subcoeff);

    multimap<double, int> norms;

    for(int i = 0; i < overlaps.mdim(); ++i) {
      const double norm = blas::dot_product(overlaps.element_ptr(0, i), overlaps.ndim(), overlaps.element_ptr(0, i));
      norms.emplace(norm, i);
    }

    cout << endl << "  o Forming dimer's active orbitals arising from " << (closed ? "closed " : "virtual ") <<  set_name
                 << " orbitals. Threshold for inclusion in cadidate space: " << setw(6) << setprecision(3) << active_thresh_ << endl;

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
      throw runtime_error("Try adjust active_thresh.");
    } else {
      const set<int> active_set(active_list.begin(), active_list.end());
      for (size_t i = 0; i < subcoeff->mdim(); ++i) {
        if (active_set.count(i)) {
          const int imo = bounds.first + i;
          assert(mask.count(imo));
          mask.erase(imo);
          copy_n(subcoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, active_position++)); //fill active columns, order=[cA,vA,cB,vB]
        }
      }
    }
  }

  //fill common closed and virtual subspace
  size_t closed_position = 0;
  for (int i = 0; i < nclosed_HF; ++i)
    if (mask.count(i))
      copy_n(treatment->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, closed_position++)); //fill closed columns

  size_t virt_position = nclosed + nact;
  for (int i = nclosed_HF; i < dimerbasis; ++i)
    if (mask.count(i))
      copy_n(treatment->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, virt_position++)); //fill virtual columns

  return out_coeff;
}


void Dimer::reduce_active(shared_ptr<const PTree> idata) {
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;
  const int nactvirtA = isolated_refs_.first->nvirt() - active_refs_.first->nvirt();
  const int nactvirtB = isolated_refs_.second->nvirt() - active_refs_.second->nvirt(); // Number of active orbitals in virtual HF subspace
  assert(nactA == isolated_refs_.first->nclosed() - active_refs_.first->nclosed() + nactvirtA);
  assert(nactB == isolated_refs_.second->nclosed() - active_refs_.second->nclosed() + nactvirtB);
  const int dimerbasis = sgeom_->nbasis();
  //sref has been activated; nclosed, nact, nvirt are defined
  const int nclosed = sref_->nclosed();
  const int nvirt_HF = sref_->nvirt() + nactvirtA + nactvirtB;
  assert(dimerbasis == nclosed + isolated_refs_.first->nclosed()-active_refs_.first->nclosed() + isolated_refs_.second->nclosed()-active_refs_.second->nclosed() + nvirt_HF);

  auto deact = idata->get_array<int,4>("reduction", {0,0,0,0});
  const int cA = deact[0];
  const int vA = deact[1];
  const int cB = deact[2];
  const int vB = deact[3];
  cout << endl << " o Active space is truncated." << endl;
  cout << "    - monomer A : " << cA << " closed and " << vA << " virtual orbitals are created from active space." << endl;
  cout << "    - monomer B : " << cB << " closed and " << vB << " virtual orbitals are created from active space." << endl;
  shared_ptr<Matrix> acoeff = sref_->coeff()->slice_copy(nclosed, nclosed+nact); // [nactA nactB]
  auto closedAB = make_shared<Matrix>(dimerbasis,cA+cB);
  auto activeAB = make_shared<Matrix>(dimerbasis,nact-cA-cB-vA-vB);
  auto virtualAB = make_shared<Matrix>(dimerbasis,vA+vB);

  //reduced to closed
  closedAB->copy_block(0,0,  dimerbasis,cA, acoeff->get_submatrix(0,0,     dimerbasis,cA));
  closedAB->copy_block(0,cA, dimerbasis,cB, acoeff->get_submatrix(0,nactA, dimerbasis,cB));
  //reduced to virtual
  virtualAB->copy_block(0,0,  dimerbasis,vA, acoeff->get_submatrix(0,nactA-vA, dimerbasis,vA));
  virtualAB->copy_block(0,vA, dimerbasis,vB, acoeff->get_submatrix(0,nact-vB,  dimerbasis,vB));
  //remains active
  activeAB->copy_block(0,0,           dimerbasis,nactA-cA-vA, acoeff->get_submatrix(0,cA,       dimerbasis,nactA-cA-vA));
  activeAB->copy_block(0,nactA-cA-vA, dimerbasis,nactB-cB-vB, acoeff->get_submatrix(0,nactA+cB, dimerbasis,nactB-cB-vB));

  shared_ptr<Matrix> outcoeff = sref_->coeff()->copy();
  outcoeff->copy_block(0,nclosed,  dimerbasis,cA+cB, closedAB);
  outcoeff->copy_block(0,nclosed+cA+cB, dimerbasis,nact-cA-cB-vA-vB, activeAB);
  outcoeff->copy_block(0,nclosed+nact-vA-vB, dimerbasis,vA+vB, virtualAB);

  //update active refs
  active_refs_ = { make_shared<Reference>(active_refs_.first->geom(), active_refs_.first->coeff(), active_refs_.first->nclosed()+cA, active_refs_.first->nact()-cA-vA, active_refs_.first->nvirt()+vA),
                   make_shared<Reference>(active_refs_.second->geom(), active_refs_.second->coeff(), active_refs_.first->nclosed()+cB, active_refs_.first->nact()-cB-vB, active_refs_.first->nvirt()+vB)};

  //synchronization
  outcoeff->broadcast();

  //semi canonical coeff
  sref_ = make_shared<Reference>(sgeom_, make_shared<Coeff>(*outcoeff), nclosed+cA+cB, nact-cA-cB-vA-vB, nvirt_HF - nactvirtA - nactvirtB + vA+vB);

}
