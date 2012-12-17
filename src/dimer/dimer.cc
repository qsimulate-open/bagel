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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <src/scf/geometry.h>
#include <src/scf/coeff.h>
#include <src/scf/matrix1e.h>
#include <src/scf/overlap.h>
#include <src/scf/fock.h>
#include <src/fci/harrison.h>
#include <src/dimer/dimer.h>

using namespace std;
using namespace bagel;

/************************************************************************************
*  Dimer::Dimer(shared_ptr<Geometry> A, vector<double> displacement)                *
*                                                                                   *
************************************************************************************/
Dimer::Dimer(shared_ptr<const Geometry> A, array<double,3> displacement) : dimerbasis_(2*A->nbasis()),
 nbasis_(A->nbasis(), A->nbasis()), symmetric_(true) {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   shared_ptr<const Geometry> geomB(new const Geometry((*A), displacement));

   geoms_ = make_pair(A, geomB);
   nele_ = make_pair(A->nele(), geomB->nele());

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();
}

Dimer::Dimer(shared_ptr<const Reference> A, array<double,3> displacement) : dimerbasis_(2*A->geom()->nbasis()),
nbasis_(A->geom()->nbasis(), A->geom()->nbasis()), symmetric_(true)
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_ = make_pair(A->coeff(), A->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB(new const Geometry((*geomA), displacement));

   geoms_ = make_pair(geomA, geomB);
   nele_ = make_pair(geomA->nele(), geomB->nele());

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();

   cout << " ===== Constructing Dimer reference ===== " << endl;
   construct_coeff(); // Constructs projected coefficients and stores them in proj_coeff;
   orthonormalize();  // Orthogonalizes projected coefficients and stores them in scoeff_;

   int nclo = 2*A->nclosed();
   int nact = 2*A->nact();
   int nvirt = 2*A->nvirt();

   shared_ptr<const Reference> tmpref(new Reference(geomB, A->coeff(), A->nclosed(), A->nact(), A->nvirt(),
            A->energy(), A->rdm1(), A->rdm2(), A->rdm1_av(), A->rdm2_av() ) );
   refs_ = make_pair(A, tmpref);

   sref_ = shared_ptr<Reference>(new Reference(sgeom_, scoeff_, nclo, nact, nvirt ));
}

Dimer::Dimer(shared_ptr<const CIWfn> A, array<double,3> displacement) : dimerbasis_(2*A->geom()->nbasis()),
nbasis_(A->geom()->nbasis(), A->geom()->nbasis()), symmetric_(true)
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_ = make_pair(A->coeff(), A->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB(new const Geometry((*geomA), displacement));

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

   shared_ptr<Reference> Aref(new Reference(geomA, coeffs_.first, ncore_.first, nact_.first, nvirt_.first));
   shared_ptr<Reference> Bref(new Reference(geomB, coeffs_.second, ncore_.second, nact_.second, nvirt_.second));
   refs_ = make_pair(Aref, Bref);

   sref_ = shared_ptr<Reference>(new Reference(sgeom_, scoeff_, nclo, nact, nvirt ));
}

void Dimer::construct_geometry() {
   vector<shared_ptr<const Geometry> > geo_vec;
   geo_vec.push_back(geoms_.first);
   geo_vec.push_back(geoms_.second);

   nbasis_ = make_pair(geoms_.first->nbasis(), geoms_.second->nbasis());

   sgeom_ = shared_ptr<Geometry>(new Geometry(geo_vec));
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

  proj_coeff_ = shared_ptr<Coeff>(new Coeff(sgeom_));
  // TODO - Ideally, these would all be projections onto the new basis.

  double *Adata = coeffs_.first->data();
  double *Bdata = coeffs_.second->data();
  double *Sdata = proj_coeff_->data();

  for(int i = 0; i < nbasisA; ++i, Adata += nbasisA) {
    Sdata = copy(Adata, Adata + nbasisA, Sdata);
    fill(Sdata, Sdata + nbasisB, 0.0); Sdata += nbasisB;
  }

  for(int i = 0; i < nbasisB; ++i, Bdata += nbasisB) {
    fill(Sdata, Sdata + nbasisA, 0.0); Sdata += nbasisA;
    Sdata = copy(Bdata, Bdata + nbasisB, Sdata);
  }

  const int ncloA = ncore_.first;
  const int ncloB = ncore_.second;

  const int nactA = nact_.first;
  const int nactB = nact_.second;

  const int nvirtA = nvirt_.first;
  const int nvirtB = nvirt_.second;

  #if 0
  shared_ptr<Matrix> tmpcoeff = proj_coeff_->slice(0,ncloA);
  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(nbasisA, nbasisA+ncloB));

  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(ncloA, ncloA+nactA));
  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(nbasisA+ncloB, nbasisA+ncloB+nactB));

  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(ncloA+nactA, ncloA+nactA+nvirtA));
  tmpcoeff = tmpcoeff->merge(proj_coeff_->slice(nbasisA+ncloB+nactB, nbasisA+ncloB+nactB+nvirtB));

  scoeff_ = shared_ptr<Coeff>(new Coeff(*tmpcoeff));
  #else
  scoeff_ = shared_ptr<Coeff>(new Coeff(*proj_coeff_));
  #endif
} 

shared_ptr<Coeff> Dimer::overlap() const {
/* What I need to do is use a sgeo to make a big overlap matrix
   and then transform it with scoeffs into the MO basis */

   /* s overlap in AO basis */
   Overlap ovlp(sgeom_);

   /* transform to MO basis with proj_coeff */
   shared_ptr<Coeff> novlp(new Coeff( (*scoeff_) % ovlp * (*scoeff_) ));

   return novlp;
}

shared_ptr<Matrix> Dimer::form_density_rhf(shared_ptr<const Coeff> coeff) const {
  shared_ptr<Matrix> out = coeff->form_density_rhf(ncore_.first);
  shared_ptr<Coeff> tmp = shared_ptr<Coeff>(new Coeff(*coeff->slice(nbasis_.first, nbasis_.first + nbasis_.second)));
  *out += *(tmp->form_density_rhf(ncore_.second));

  return out;
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

   scoeff_ = shared_ptr<Coeff>(new Coeff(*scoeff_ * S_1_2));
}

void Dimer::fci(multimap<string,string> idata) {
  { // Start FCI on unit A
    // Move occupied orbitals of unit B to form the core orbitals
    shared_ptr<Matrix> Amatrix = scoeff_->slice(nbasis_.first, nbasis_.first + ncore_.second)->merge(scoeff_->slice(0,nbasis_.first));
    shared_ptr<Coeff> Acoeff(new Coeff(*Amatrix));

    // Set up variables for this fci
    shared_ptr<Reference> Aref(new Reference(sgeom_, Acoeff, (nele_.first + nele_.second)/2, 0, 0));

    const int ncore = ncore_.second; // ncore in this set up is all the occupied orbitals of unit B
    const int norb  = nbasis_.first; // Assuming use of the whole space

    // for now let's just hardcode in HZ. This can maybe be changed later
    shared_ptr<FCI> fci(new HarrisonZarrabian(idata, Aref, ncore, norb));

    fci->compute();
    ccvecs_.first = fci->civectors();
  }

  { // Start FCI on unit B
    shared_ptr<Matrix> Bmatrix = scoeff_->slice(0,ncore_.first)->merge(scoeff_->slice(nbasis_.first, nbasis_.first + nbasis_.second));
    shared_ptr<Coeff> Bcoeff(new Coeff(*Bmatrix));

    // Set up variables for this fci
    shared_ptr<Reference> Bref(new Reference(sgeom_, Bcoeff, (nele_.first + nele_.second)/2, 0, 0));

    const int ncore = ncore_.first;
    const int norb = nbasis_.second;

    // HZ fci
    shared_ptr<FCI> fci(new HarrisonZarrabian(idata, Bref, ncore, norb));

    fci->compute();
    ccvecs_.second = fci->civectors();
  }

  // Wavefunctions came from FCI calculations, so I should know the following
  ncore_ = make_pair(0,0);
  nact_ = make_pair(nbasis_.first, nbasis_.second);
  nvirt_ = make_pair(0,0);

  symmetric_ = false;
}

void Dimer::localize(multimap<string, string> idata) {
  string localizemethod = read_input<string>(idata,"localization", "pm");

  shared_ptr<OrbitalLocalization> localization;
  if (localizemethod == "region") {
    vector<int> sizes = { geoms_.first->natom(), geoms_.second->natom() };
    localization = shared_ptr<OrbitalLocalization>(new RegionLocalization(sref_, sizes));
  }
  else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey") {
    localization = shared_ptr<OrbitalLocalization>(new PMLocalization(sref_));
  }
  else throw std::runtime_error("Unrecognized orbital localization method");

  set_coeff(localization->localize());
}
