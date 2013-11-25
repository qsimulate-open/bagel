//
// BAGEL - Parallel electron correlation program.
// Filename: dimer.cc
// Copyright (C) 2012 Shane Parker
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
#include <src/wfn/geometry.h>
#include <src/molecule/overlap.h>
#include <src/scf/coeff.h>
#include <src/scf/fock.h>
#include <src/fci/harrison.h>
#include <src/fci/knowles.h>
#include <src/wfn/reference.h>
#include <src/molecule/localization.h>
#include <src/util/lexical_cast.h>

using namespace std;
using namespace bagel;

/************************************************************************************
*  Dimer::Dimer(shared_ptr<Geometry> A, vector<double> displacement)                *
*                                                                                   *
************************************************************************************/
Dimer::Dimer(shared_ptr<const Geometry> A, array<double,3> displacement) : dimerbasis_(2*A->nbasis()),
 nbasis_(A->nbasis(), A->nbasis()) {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   auto geomB = make_shared<const Geometry>((*A), displacement);

   geoms_ = make_pair(A, geomB);
   nele_ = make_pair(A->nele(), geomB->nele());

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();
}

Dimer::Dimer(shared_ptr<const Reference> A, array<double,3> displacement) : dimerbasis_(2*A->geom()->nbasis()),
nbasis_(A->geom()->nbasis(), A->geom()->nbasis())
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_ = make_pair(A->coeff(), A->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB = make_shared<const Geometry>((*geomA), displacement);

   geoms_ = make_pair(geomA, geomB);
   nele_ = make_pair(geomA->nele(), geomB->nele());

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();

   cout << " ===== Constructing Dimer reference ===== " << endl;
   construct_coeff(); // Constructs projected coefficients and stores them in proj_coeff;
   orthonormalize();  // Orthogonalizes projected coefficients and stores them in scoeff_;

   nclosed_ = 2*A->nclosed();
   int nact = 2*A->nact();
   int nvirt = 2*A->nvirt();

   auto tmpref = make_shared<const Reference>(geomB, A->coeff(), A->nclosed(), A->nact(), A->nvirt(),
            A->energy(), A->rdm1(), A->rdm2(), A->rdm1_av(), A->rdm2_av() );
   refs_ = make_pair(A, tmpref);

   sref_ = make_shared<Reference>(sgeom_, scoeff_, nclosed_, nact, nvirt);
}

Dimer::Dimer(shared_ptr<const Reference> A, shared_ptr<const Reference> B) : dimerbasis_(A->geom()->nbasis() + B->geom()->nbasis()),
nbasis_(A->geom()->nbasis(), B->geom()->nbasis())
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_ = make_pair(A->coeff(), B->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB = B->geom();

   geoms_ = make_pair(geomA, geomB);
   nele_ = make_pair(geomA->nele(), geomB->nele());

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();

   cout << " ===== Constructing Dimer reference ===== " << endl;
   construct_coeff(); // Constructs projected coefficients and stores them in proj_coeff;
   orthonormalize();  // Orthogonalizes projected coefficients and stores them in scoeff_;

   nclosed_ = A->nclosed() + B->nclosed();
   int nact = A->nact() + B->nact();
   int nvirt = A->nvirt() + B->nact();

   refs_ = make_pair(A, B);

   sref_ = make_shared<Reference>(sgeom_, scoeff_, nclosed_, nact, nvirt);
}

Dimer::Dimer(shared_ptr<const Reference> superref, pair<int,int> regions) : sgeom_(superref->geom()), dimerbasis_(superref->geom()->nbasis())
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
  int nbasisA = 0;
  int neleA = 0;
  for (int i = 0; i < regions.first; ++i) {
    const auto& atom = sgeom_->atoms(i);
    nbasisA += atom->nbasis();
    neleA += atom->atom_number();
  }

  int nbasisB = 0;
  int neleB = 0;
  for (int i = regions.first; i < regions.first + regions.second; ++i) {
    const auto& atom = sgeom_->atoms(i);
    nbasisB += atom->nbasis();
    neleB += atom->atom_number();
  }

  assert(nbasisA + nbasisB == dimerbasis_);
  nbasis_ = make_pair(nbasisA, nbasisB);

  const int nele = neleA + neleB;
  nele_ = make_pair(neleA, neleB);

  const int nocc = nele/2;
  const int noccA = neleA/2;
  const int noccB = neleB/2;

  const int nclosed = superref->nclosed();
  const int nact = superref->nact();

  vector<int> sizes = {{regions.first, regions.second}}; // a little bit of a hack but can be improved later
  shared_ptr<Matrix> tmpcoeff = superref->coeff()->slice(nclosed, nclosed + nact);
  auto localization = make_shared<RegionLocalization>(make_shared<PTree>(), sgeom_, tmpcoeff, sizes, nact, 0, 0);
  auto local_active = localization->localize();
  auto local_coeff = make_shared<Coeff>(*superref->coeff());
  copy_n(local_active->element_ptr(0, 0), dimerbasis_*nact, local_coeff->element_ptr(0,nclosed));

  nclosed_ = nclosed;
  vector<int> act_regions = localization->region_orbitals(0);
  nact_ = make_pair(act_regions.at(0), act_regions.at(1));
  //nfilledactive_ = make_pair(noccA - ncore_.first, noccB - ncore_.second);

  auto density = superref->coeff()->form_density_rhf(nocc);

  shared_ptr<const Matrix> hcore(new Hcore(sgeom_));
  shared_ptr<Matrix> fock(new Fock<1>(sgeom_, hcore, density, superref->coeff()));
  Matrix intermediate((*local_coeff) % (*fock) * (*local_coeff));

  Matrix transform(dimerbasis_, dimerbasis_); transform.unit();

  const int subsize = nclosed_ + nact_.first + nact_.second;
  vector<int> subsizes = {nclosed_, nact_.first, nact_.second};
  vector<double> subeigs(dimerbasis_, 0.0);

  shared_ptr<Matrix> diag_blocks = intermediate.get_submatrix(0, 0, subsize, subsize)->diagonalize_blocks(subeigs.data(), subsizes);
  transform.copy_block(0,0,subsize,subsize,diag_blocks);

  multimap<double, int> active_fock;
  for (int i = nclosed_; i < nact + nclosed_; ++i) {
    active_fock.insert(make_pair(subeigs.at(i), i));
  }

  int filledA = 0;
  int filledB = 0;
  auto iactive = active_fock.begin();
  for ( int totalfilledactive = nocc - nclosed_; totalfilledactive != 0; --totalfilledactive, ++iactive ) {
    int imo = iactive->second;
    if ( (nclosed_ <= imo) && (imo < nclosed_ + nact_.first) ) ++filledA;
    else ++filledB;
  }
  nfilledactive_ = make_pair(filledA, filledB);

  shared_ptr<Matrix> ncoeff(new Matrix(*local_coeff * transform));
  scoeff_ = make_shared<Coeff>(*ncoeff);
  sref_ = make_shared<Reference>(superref, scoeff_);
  sref_->set_eig(subeigs);
}

Dimer::Dimer(shared_ptr<const CIWfn> A, array<double,3> displacement) : dimerbasis_(2*A->geom()->nbasis()),
nbasis_(A->geom()->nbasis(), A->geom()->nbasis())
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_ = make_pair(A->coeff(), A->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB = make_shared<const Geometry>((*geomA), displacement);

   geoms_ = make_pair(geomA, geomB);
   nele_ = make_pair(geomA->nele(), geomB->nele());

   ci_ = make_pair(A, A);
   ccvecs_ = make_pair(A->civectors(), A->civectors());

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();

   cout << " ===== Constructing Dimer reference ===== " << endl;
   construct_coeff(); // Constructs projected coefficients and stores them in proj_coeff;
   orthonormalize();  // Orthogonalizes projected coefficients and stores them in scoeff_;

   int nclo = 2*A->ncore();
   int nact = 2*A->nact();
   int nvirt = 2*A->nvirt();

   auto Aref = make_shared<Reference>(geomA, coeffs_.first, ncore_.first, nact_.first, nvirt_.first);
   auto Bref = make_shared<Reference>(geomB, coeffs_.second, ncore_.second, nact_.second, nvirt_.second);
   refs_ = make_pair(Aref, Bref);

   sref_ = make_shared<Reference>(sgeom_, scoeff_, nclo, nact, nvirt);
}

void Dimer::construct_geometry() {
   vector<shared_ptr<const Geometry>> geo_vec;
   geo_vec.push_back(geoms_.first);
   geo_vec.push_back(geoms_.second);

   nbasis_ = make_pair(geoms_.first->nbasis(), geoms_.second->nbasis());

   sgeom_ = make_shared<Geometry>(geo_vec);
}

void Dimer::construct_coeff() {
  const int nbasisA = nbasis_.first;
  const int nbasisB = nbasis_.second;

  if(static_cast<bool>(refs_.first)) {
    ncore_.first  = refs_.first->nclosed();
    ncore_.second = refs_.second->nclosed();

    nact_.first  = refs_.first->nact();
    nact_.second = refs_.second->nact();

    nvirt_.first  = refs_.first->nvirt();
    nvirt_.second = refs_.second->nvirt();
  }
  else if (static_cast<bool>(ci_.first)) {
    ncore_.first  = ci_.first->ncore();
    ncore_.second = ci_.second->ncore();

    nact_.first  = ci_.first->nact();
    nact_.second = ci_.second->nact();

    nvirt_.first  = ci_.first->nvirt();
    nvirt_.second = ci_.second->nvirt();
  }
  else {
    // Round nele up for number of orbitals
    ncore_.first  = (nele_.first + 1)/2;
    ncore_.second = (nele_.second + 1)/2;

    nact_.first  = 0;
    nact_.second = 0;

    nvirt_.first = (nbasisA - ncore_.first);
    nvirt_.second = (nbasisB - ncore_.second);
  }

  proj_coeff_ = make_shared<Coeff>(sgeom_);
  // TODO - Ideally, these would all be projections onto the new basis.

  const double *Adata = coeffs_.first->data();
  const double *Bdata = coeffs_.second->data();
  double *Sdata = proj_coeff_->data();

  for(int i = 0; i < nbasisA; ++i, Adata += nbasisA) {
    Sdata = copy_n(Adata, nbasisA, Sdata);
    fill_n(Sdata, nbasisB, 0.0); Sdata += nbasisB;
  }

  for(int i = 0; i < nbasisB; ++i, Bdata += nbasisB) {
    fill_n(Sdata, nbasisA, 0.0); Sdata += nbasisA;
    Sdata = copy_n(Bdata, nbasisB, Sdata);
  }

  const int ncloA = ncore_.first;
  const int ncloB = ncore_.second;

  const int nactA = nact_.first;
  const int nactB = nact_.second;

  const int nvirtA = nvirt_.first;
  const int nvirtB = nvirt_.second;

  shared_ptr<Matrix> tmpcoeff = proj_coeff_->slice(0,ncloA);
  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(nbasisA, nbasisA+ncloB));

  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(ncloA, ncloA+nactA));
  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(nbasisA+ncloB, nbasisA+ncloB+nactB));

  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(ncloA+nactA, ncloA+nactA+nvirtA));
  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(nbasisA+ncloB+nactB, nbasisA+ncloB+nactB+nvirtB));

  scoeff_ = make_shared<Coeff>(*tmpcoeff);
}

shared_ptr<Coeff> Dimer::overlap() const {
   Overlap ovlp(sgeom_);
   return make_shared<Coeff>((*scoeff_) % ovlp * (*scoeff_));
}

void Dimer::orthonormalize() {
   shared_ptr<Coeff> ovlp = overlap();
   Matrix S = *ovlp;

   unique_ptr<double[]> eig(new double[dimerbasis_]);
   S.diagonalize(eig.get());

   Matrix S_1_2(S);
   double *S12_data = S_1_2.data();
   for( int ii = 0; ii < dimerbasis_; ++ii) {
      dscal_(dimerbasis_, 1.0/sqrt(eig[ii]), S12_data, 1);
      S12_data += dimerbasis_;
   }

   S_1_2 = S_1_2 ^ S;

   scoeff_ = make_shared<Coeff>(*scoeff_ * S_1_2);
}

void Dimer::embed_refs() {
  const int noccA = nele_.first/2;
  const int noccB = nele_.second/2;
  const int nocc  = noccA + noccB;

  const int nclosed = nclosed_;

  // filled_active is the number of orbitals in the active space that should be filled
  const int filled_activeA = nfilledactive_.first;
  const int filled_activeB = nfilledactive_.second;

  const int nactA = nact_.first;
  const int nactB = nact_.second;
  const int nact = nactA + nactB;

  const int nbasisA = nbasis_.first;
  const int nbasisB = nbasis_.second;

  { // Move occupied orbitals of unit B to form the core orbitals
    auto Amatrix = make_shared<Matrix>(dimerbasis_, dimerbasis_);
    Amatrix->copy_block(0, 0, dimerbasis_, nclosed, scoeff_->element_ptr(0,0)); // Total closed space
    Amatrix->copy_block(0, nclosed, dimerbasis_, filled_activeB, scoeff_->element_ptr(0,nclosed + nactA)); // FilledActive B
    Amatrix->copy_block(0, nclosed + filled_activeB, dimerbasis_, nactA, scoeff_->element_ptr(0,nclosed)); // Active A
    auto Acoeff = make_shared<Coeff>(*Amatrix);

    // Set up variables for this fci
    const int ncore = nclosed + filled_activeB;
    const int norb  = nactA;

    embedded_refs_.first = make_shared<Reference>(sgeom_, Acoeff, ncore, norb, 0);
  }

  { // Move occupied orbitals of unit A to form core of unit B
    auto Bmatrix = make_shared<Matrix>(dimerbasis_, dimerbasis_);
    Bmatrix->copy_block(0, 0, dimerbasis_, nclosed, scoeff_->element_ptr(0,0)); // Total closed space
    Bmatrix->copy_block(0, nclosed, dimerbasis_, filled_activeA, scoeff_->element_ptr(0,nclosed)); // FilledActive A
    Bmatrix->copy_block(0, nclosed + filled_activeA, dimerbasis_, nactB, scoeff_->element_ptr(0,nclosed + nactA)); // Active B
    auto Bcoeff = make_shared<Coeff>(*Bmatrix);

    // Set up variables for this fci
    const int ncore = nclosed + filled_activeA;
    const int norb  = nactB;

    embedded_refs_.second = make_shared<Reference>(sgeom_, Bcoeff, ncore, norb, 0);
  }
}

void Dimer::localize(const std::shared_ptr<const PTree> idata) {
  string localizemethod = idata->get<string>("algorithm", "pm");

  shared_ptr<OrbitalLocalization> localization;
  if (localizemethod == "region") {
    vector<int> sizes = { geoms_.first->natom(), geoms_.second->natom() };
    localization = make_shared<RegionLocalization>(idata, sref_, sizes);
  }
  else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey") {
    localization = make_shared<PMLocalization>(idata, sref_);
  }
  else throw std::runtime_error("Unrecognized orbital localization method");

  shared_ptr<const Matrix> local_coeff = localization->localize();
  auto S = make_shared<Overlap>(sgeom_);
  auto overlaps = make_shared<Matrix>((*proj_coeff_) % (*S) * (*local_coeff));

  const int nclosed = nclosed_;
  const int nact = nact_.first + nact_.second;
  const int nvirt = nvirt_.first + nvirt_.second;

  multimap<double, int> Anorms;
  set<int> closed_setA, closed_setB;
  for(int i = 0; i < nclosed; ++i) {
    double* cdata = overlaps->element_ptr(0,i);
    double norm = ddot_(nbasis_.first, cdata, 1, cdata, 1);

    if (norm > 0.7) closed_setA.insert(i);
    else if (norm < 0.3) closed_setB.insert(i);
    else throw runtime_error("Trouble in classifying orbitals");
  }

  set<int> active_setA, active_setB;
  for(int i = nclosed; i < nclosed + nact; ++i) {
    double* cdata = overlaps->element_ptr(0,i);
    double norm = ddot_(nbasis_.first, cdata, 1, cdata, 1);

    if (norm > 0.7) active_setA.insert(i);
    else if (norm < 0.3) active_setB.insert(i);
    else throw runtime_error("Trouble in classifying orbitals");
  }

  auto new_coeff = make_shared<Matrix>(dimerbasis_, dimerbasis_);
  int imo = 0;

  for(auto& iset : closed_setA) {
    copy_n(local_coeff->element_ptr(0,iset), dimerbasis_, new_coeff->element_ptr(0,imo));
    ++imo;
  }
  for(auto& iset : closed_setB) {
    copy_n(local_coeff->element_ptr(0,iset), dimerbasis_, new_coeff->element_ptr(0,imo));
    ++imo;
  }

  for(auto& iset : active_setA) {
    copy_n(local_coeff->element_ptr(0,iset), dimerbasis_, new_coeff->element_ptr(0,imo));
    ++imo;
  }
  for(auto& iset : active_setB) {
    copy_n(local_coeff->element_ptr(0,iset), dimerbasis_, new_coeff->element_ptr(0,imo));
    ++imo;
  }

  copy_n(local_coeff->element_ptr(0,nclosed + nact), dimerbasis_*nvirt, new_coeff->element_ptr(0,imo));

  auto out = make_shared<const Coeff>(*new_coeff);

  set_coeff(out);
}

void Dimer::set_active(const std::shared_ptr<const PTree> idata) {
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
  pair<shared_ptr<const Reference>, shared_ptr<const Reference>> active_refs =
        make_pair(refs_.first->set_active(Alist), refs_.second->set_active(Blist));

  // Update Dimer info
  const int nclosedA = active_refs.first->nclosed();
  const int nclosedB = active_refs.second->nclosed();
  nclosed_ = nclosedA + nclosedB;

  const int nactA = active_refs.first->nact();
  const int nactB = active_refs.second->nact();
  nact_ = make_pair(nactA, nactB);
  const int nact = nactA + nactB;

  const int nvirtA = active_refs.first->nvirt();
  const int nvirtB = active_refs.second->nvirt();
  nvirt_ = make_pair(nvirtA, nvirtB);

  auto active = make_shared<Matrix>(dimerbasis_, nact);

  const int nbasisA = nbasis_.first;
  const int nbasisB = nbasis_.second;

  active->copy_block(0, 0, nbasisA, nactA, active_refs.first->coeff()->get_block(0, nclosedA, nbasisA, nactA));
  active->copy_block(nbasisA, nactA, nbasisB, nactB, active_refs.second->coeff()->get_block(0, nclosedB, nbasisB, nactB));

  auto S = make_shared<Overlap>(sgeom_);
  auto overlaps = make_shared<Matrix>((*sref_->coeff()) % (*S) * (*active));
  double* odata = overlaps->data();

  multimap<double, int> norms;

  for(int i = 0; i < dimerbasis_; ++i, ++odata) {
    double norm = ddot_(nact, odata, dimerbasis_, odata, dimerbasis_);
    norms.insert(make_pair(norm, i));
  }

  auto norm_iter = norms.rbegin();
  set<int> active_list;
  for(int i = 0; i < nact; ++i, ++norm_iter) active_list.insert(norm_iter->second);

  auto out = sref_->set_active(active_list);
  const int nfilledA = geoms_.first->nele()/2 - nclosedA;
  const int nfilledB = geoms_.second->nele()/2 - nclosedB;
  nfilledactive_ = make_pair( nfilledA, nfilledB );

  set_sref(out);
}

// RHF and then localize
void Dimer::scf(const shared_ptr<const PTree> idata) {
  Timer dimertime;

  // SCF
  auto hfdata = idata->get_child_optional("hf") ? idata->get_child_optional("hf") : make_shared<const PTree>();
  shared_ptr<SCF> rhf = dynamic_pointer_cast<SCF>(construct_method("hf", hfdata, sgeom_, sref_));
  rhf->compute();
  set_sref(rhf->conv_to_ref());
  dimertime.tick_print("Dimer SCF");

  shared_ptr<Matrix> dimerdensity = sref_->coeff()->form_density_rhf(nclosed_);
  shared_ptr<Matrix> dimercoeff = scoeff_->slice(0,nclosed_);

  // Set active space based on overlap
  if (proj_coeff_ == nullptr) throw runtime_error("For Dimer::driver, Dimer must be constructed from a HF reference");
  else set_active(idata);

  // Localize
  const string localmethod = idata->get<string>("localization", "default");
  dimertime.tick();
  if (localmethod != "none") {
    shared_ptr<const PTree> localize_data = idata->get_child_optional("localization");
    if (!localize_data) localize_data = make_shared<const PTree>();

    localize(localize_data);
    dimertime.tick_print("Dimer localization");

    // Sub-diagonalize Fock Matrix
    auto fock  = make_shared<const Fock<1>>(sgeom_, sref_->hcore(), dimerdensity, dimercoeff);
    Matrix intermediate((*scoeff_) % (*fock) * (*scoeff_));
    dimertime.tick_print("Dimer Fock matrix formation");

    Matrix transform(dimerbasis_, dimerbasis_); transform.unit();

    const int subsize = nclosed_ + nact_.first + nact_.second;
    vector<int> subsizes = {nclosed_, nact_.first, nact_.second};
    vector<double> subeigs(dimerbasis_, 0.0);

    shared_ptr<Matrix> diag_blocks = intermediate.get_submatrix(0, 0, subsize, subsize)->diagonalize_blocks(subeigs.data(), subsizes);
    transform.copy_block(0,0,subsize,subsize,diag_blocks);

    auto ncoeff = make_shared<Matrix>(*scoeff_ * transform);
    set_coeff(ncoeff);
    sref_->set_eig(subeigs);
    dimertime.tick_print("Fock block diagonalization");
  }
}


shared_ptr<DimerCAS> Dimer::compute_cispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(nfilledactive().first, nfilledactive().second);
  pair<int,int> neleb = make_pair(nfilledactive().first, nfilledactive().second);

  auto d1 = make_shared<Determinants>(nact().first, nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = make_shared<Determinants>(nact().second, nelea.first, neleb.first, /*compress*/false, /*mute*/true);
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
  pair<int,int> nelea = make_pair(nfilledactive().first, nfilledactive().second);
  pair<int,int> neleb = make_pair(nfilledactive().first, nfilledactive().second);

  auto d1 = make_shared<Determinants>(nact().first, nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = make_shared<Determinants>(nact().second, nelea.first, neleb.first, /*compress*/false, /*mute*/true);
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
  pair<int,int> nelea = make_pair(nfilledactive().first, nfilledactive().second);
  pair<int,int> neleb = make_pair(nfilledactive().first, nfilledactive().second);

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
