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
#include <src/molecule/overlap.h>

using namespace std;
using namespace bagel;

/************************************************************************************
 *  Single reference plus translation vector constructors                            *
 ************************************************************************************/
Dimer::Dimer(shared_ptr<const PTree> input, shared_ptr<const Geometry> A) : input_(input) {
    array<double, 3> translation = input->get_array<double, 3>("translate");
    if (input->get<bool>("angstrom", false))
      for_each(translation.begin(), translation.end(), [] (double& p) { p/= au2angstrom__; });
    auto geomB = make_shared<const Geometry>((*A), translation);

    geoms_ = make_pair(A, geomB);
    construct_geometry();
  }

Dimer::Dimer(shared_ptr<const PTree> input, shared_ptr<const Reference> A) : input_(input) {
  array<double, 3> translation = input->get_array<double, 3>("translate");
  if (input->get<bool>("angstrom", false))
    for_each(translation.begin(), translation.end(), [] (double& p) { p/= au2angstrom__; });

  assert(A);
  auto geomB = make_shared<const Geometry>((*A->geom()), translation);
  geoms_ = make_pair(A->geom(), geomB);
  construct_geometry();

  auto tmpref = make_shared<const Reference>(geomB, A->coeff(), A->nclosed(), A->nact(), A->nvirt(),
      A->energy(), A->rdm1(), A->rdm2(), A->rdm1_av(), A->rdm2_av() );
  isolated_refs_ = make_pair(A, tmpref);
  shared_ptr<const Matrix> coeff = construct_coeff();

  sref_ = make_shared<Reference>(sgeom_, make_shared<const Coeff>(move(*coeff)), 2*A->nclosed(), 2*A->nact(), 2*A->nvirt());
}

Dimer::Dimer(shared_ptr<const PTree> input, shared_ptr<const Reference> A, shared_ptr<const Reference> B) : input_(input) {
  geoms_ = make_pair(A->geom(), B->geom());
  construct_geometry();

  isolated_refs_ = make_pair(A, B);
  shared_ptr<const Matrix> coeff = construct_coeff();

  sref_ = make_shared<Reference>(sgeom_, make_shared<const Coeff>(move(*coeff)), A->nclosed() + B->nclosed(), A->nact() + B->nact(), A->nvirt() + B->nvirt());
}


void Dimer::construct_geometry() {
  cout << " ===== Constructing Dimer geometry ===== " << endl;

  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  if (mdata) {
    Muffle hide_cout;
    geoms_ = make_pair(make_shared<Geometry>(*geoms_.first, mdata), make_shared<Geometry>(*geoms_.second, mdata));
  }

  vector<shared_ptr<const Geometry>> geo_vec = {{ geoms_.first, geoms_.second }};
  shared_ptr<const PTree> env_data = input_->get_child_optional("environment");
  if (env_data) {
    Muffle hide_cout;
    geo_vec.push_back(make_shared<Geometry>(env_data));
  }
  sgeom_ = make_shared<Geometry>(geo_vec);
}

shared_ptr<const Matrix> Dimer::form_projected_coeffs() {
  shared_ptr<const Matrix> projectedA = isolated_refs_.first->project_coeff(sgeom_)->coeff();
  shared_ptr<const Matrix> projectedB = isolated_refs_.second->project_coeff(sgeom_)->coeff();

  shared_ptr<Matrix> projected = projectedA->merge(projectedB);

  const int ncloA = isolated_refs_.first->nclosed();
  const int ncloB = isolated_refs_.second->nclosed();

  const int nactA = isolated_refs_.first->nact();
  const int nactB = isolated_refs_.second->nact();

  const int nvirtA = isolated_refs_.first->nvirt();
  const int nvirtB = isolated_refs_.second->nvirt();

  assert(isolated_refs_.first->coeff()->mdim()  == ncloA + nactA + nvirtA);
  assert(isolated_refs_.second->coeff()->mdim() == ncloB + nactB + nvirtB);

  const size_t Amos = isolated_refs_.first->coeff()->mdim();

  // form "projected" coefficients
  const int dimerbasis = sgeom_->nbasis();
  auto out = projected->clone();

  size_t current = 0;
  auto cp_block = [&current, &out, &dimerbasis] (const size_t msize, const double* source) {
    out->copy_block(0, current, dimerbasis, msize, source); current += msize;
  };

  cp_block(ncloA, projected->element_ptr(0,0));
  cp_block(ncloB, projected->element_ptr(0, Amos));
  cp_block(nactA, projected->element_ptr(0, ncloA));
  cp_block(nactB, projected->element_ptr(0, Amos + ncloB));
  cp_block(nvirtA, projected->element_ptr(0, ncloA + nactA));
  cp_block(nvirtB, projected->element_ptr(0, Amos + ncloB + nactB));

  return out;
}

shared_ptr<const Matrix> Dimer::construct_coeff() {
  cout << " ===== Constructing Dimer reference =====" << endl;

  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  if (mdata) {
    isolated_refs_ = make_pair(isolated_refs_.first->project_coeff(geoms_.first), isolated_refs_.second->project_coeff(geoms_.second));
  }

  shared_ptr<const Matrix> projected = form_projected_coeffs();

  // orthonormalize the "projected" coefficients
  shared_ptr<Matrix> S = make_shared<Overlap>(sgeom_);
  auto S_invhalf = make_shared<Matrix>((*projected) % (*S) * (*projected));
  S_invhalf->inverse_half();

  return make_shared<Matrix>(*projected * *S_invhalf);
}

void Dimer::embed_refs() {
  Timer timer;
  const int nclosed = sref_->nclosed();

  // filled_active is the number of orbitals in the active space that should be filled
  const int filled_activeA = isolated_refs_.first->nclosed() - active_refs_.first->nclosed();
  const int filled_activeB = isolated_refs_.second->nclosed() - active_refs_.second->nclosed();

  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();

  shared_ptr<const Matrix> scoeff = sref_->coeff();

  const int dimerbasis = sgeom_->nbasis();

  { // Move occupied orbitals of unit B to form the core orbitals
    auto Amatrix = make_shared<Matrix>(dimerbasis, nclosed + filled_activeB + nactA);
    Amatrix->copy_block(0, 0, dimerbasis, nclosed, scoeff->element_ptr(0,0)); // Total closed space
    Amatrix->copy_block(0, nclosed, dimerbasis, filled_activeB, scoeff->element_ptr(0,nclosed + nactA)); // FilledActive B
    Amatrix->copy_block(0, nclosed + filled_activeB, dimerbasis, nactA, scoeff->element_ptr(0,nclosed)); // Active A
    auto Acoeff = make_shared<Coeff>(move(*Amatrix));

    // Set up variables for this fci
    const int ncore = nclosed + filled_activeB;
    const int norb  = nactA;

    embedded_refs_.first = make_shared<Reference>(sgeom_, Acoeff, ncore, norb, 0);
  }

  { // Move occupied orbitals of unit A to form core of unit B
    auto Bmatrix = make_shared<Matrix>(dimerbasis, nclosed + filled_activeA + nactB);
    Bmatrix->copy_block(0, 0, dimerbasis, nclosed, scoeff->element_ptr(0,0)); // Total closed space
    Bmatrix->copy_block(0, nclosed, dimerbasis, filled_activeA, scoeff->element_ptr(0,nclosed)); // FilledActive A
    Bmatrix->copy_block(0, nclosed + filled_activeA, dimerbasis, nactB, scoeff->element_ptr(0,nclosed + nactA)); // Active B
    auto Bcoeff = make_shared<Coeff>(move(*Bmatrix));

    // Set up variables for this fci
    const int ncore = nclosed + filled_activeA;
    const int norb  = nactB;

    embedded_refs_.second = make_shared<Reference>(sgeom_, Bcoeff, ncore, norb, 0);
  }
}

void Dimer::get_spaces(shared_ptr<const PTree> idata, vector<vector<int>>& spaces_A, vector<vector<int>>& spaces_B) {
  auto space = idata->get_child_optional("space");
  if (space) {
    for (auto& s : *space) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    for (auto& s : *spacea) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    for (auto& s : *spaceb) { spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
  }
}
