//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimer.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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
#include <src/mat1e/overlap.h>

using namespace std;
using namespace bagel;

/************************************************************************************
 *  Single reference plus translation vector constructors                            *
 ************************************************************************************/
Dimer::Dimer(shared_ptr<const PTree> input, shared_ptr<const Geometry> A) : input_(input) {
  array<double, 3> translation = input->get_array<double, 3>("translate");
  if (input->get<bool>("angstrom", false))
    for_each(translation.begin(), translation.end(), [](double& p) { p/= au2angstrom__; });
  auto geomB = make_shared<const Geometry>((*A), translation);

  geoms_ = {A, geomB};
  construct_geometry();
}

Dimer::Dimer(shared_ptr<const PTree> input, shared_ptr<const Reference> A, const bool linked) : input_(input) {
  if (!linked) {
    array<double, 3> translation = input->get_array<double, 3>("translate");
    if (input->get<bool>("angstrom", false))
      for_each(translation.begin(), translation.end(), [](double& p) { p/= au2angstrom__; });

    assert(A);
    auto geomB = make_shared<const Geometry>((*A->geom()), translation);
    geoms_ = {A->geom(), geomB};
    construct_geometry();

    auto tmpref = make_shared<const Reference>(geomB, A->coeff(), A->nclosed(), A->nact(), A->nvirt(),
                                               A->energy(), A->rdm1(), A->rdm2(), A->rdm1_av(), A->rdm2_av());
    isolated_refs_ = {A, tmpref};
    shared_ptr<const Matrix> coeff = construct_coeff();

    nvirt_ = {A->nvirt(), A->nvirt()};
    sref_ = make_shared<Reference>(sgeom_, make_shared<const Coeff>(move(*coeff)), 2*A->nclosed(), 2*A->nact(), 2*A->nvirt());
  }
  else { //covalently linked dimer
    print_orbital_ = input->get<bool>("print_orbital", false);
    geoms_ = {A->geom(), A->geom()}; //initialize
    construct_geometry(linked); //geoms_ is redefined with new basis sets, sgeom_ is defined
    isolated_refs_ = {A, A}; //initialize
    shared_ptr<const Matrix> coeff = construct_coeff(linked); //isolated_refs_ is redifined with new basis sets, dimension = (large,small)

    nvirt_ = {A->nvirt(), A->nvirt()};
    sref_ = make_shared<Reference>(sgeom_, make_shared<const Coeff>(move(*coeff)), A->nclosed(), /*0*/A->nact(), A->nvirt());
  }
}

Dimer::Dimer(shared_ptr<const PTree> input, shared_ptr<const Reference> A, shared_ptr<const Reference> B) : input_(input) {
  geoms_ = {A->geom(), B->geom()};
  construct_geometry();

  isolated_refs_ = {A, B};
  shared_ptr<const Matrix> coeff = construct_coeff();

  nvirt_ = {A->nvirt(), B->nvirt()};
  sref_ = make_shared<Reference>(sgeom_, make_shared<const Coeff>(move(*coeff)), A->nclosed()+B->nclosed(), A->nact()+B->nact(), A->nvirt()+B->nvirt());
}


void Dimer::construct_geometry(bool linked) {
  cout << " ===== Constructing Dimer geometry ===== " << endl;

  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  if (mdata) {
    Muffle hide_cout;
    geoms_ = {make_shared<Geometry>(*geoms_.first, mdata), make_shared<Geometry>(*geoms_.second, mdata)};
  }

  vector<shared_ptr<const Geometry>> geo_vec;
  if (!linked)
    geo_vec = { geoms_.first, geoms_.second };
  else
    geo_vec = { geoms_.first }; //linked dimer geometry, only one is required

  shared_ptr<const PTree> env_data = input_->get_child_optional("environment");
  if (env_data) {
    Muffle hide_cout;
    geo_vec.push_back(make_shared<Geometry>(env_data));
  }
  sgeom_ = make_shared<Geometry>(geo_vec);
}

shared_ptr<const Matrix> Dimer::form_projected_coeffs() {
  shared_ptr<const Matrix> projectedA = isolated_refs_.first->project_coeff(sgeom_, false)->coeff();
  shared_ptr<const Matrix> projectedB = isolated_refs_.second->project_coeff(sgeom_, false)->coeff();

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
  auto out = projected->clone();

  size_t current = 0;
  auto cp_block = [&current, &out](const MatView& source) {
    out->copy_block(0, current, source.ndim(), source.mdim(), source);
    current += source.mdim();
  };

  cp_block(projected->slice(               0,      ncloA));
  cp_block(projected->slice(Amos            , Amos+ncloB));
  cp_block(projected->slice(     ncloA      ,      ncloA+nactA));
  cp_block(projected->slice(Amos+ncloB      , Amos+ncloB+nactB));
  cp_block(projected->slice(     ncloA+nactA,      ncloA+nactA+nvirtA));
  cp_block(projected->slice(Amos+ncloB+nactB, Amos+ncloB+nactB+nvirtB));

  return out;
}


shared_ptr<const Matrix> Dimer::construct_coeff(const bool linked) {
  cout << " ===== Constructing Dimer reference =====" << endl;

  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  if (mdata) {
    isolated_refs_ = {isolated_refs_.first->project_coeff(geoms_.first, false), isolated_refs_.second->project_coeff(geoms_.second, false)};
  }

  shared_ptr<const Matrix> projected;
  if (!linked)
    projected = form_projected_coeffs();
  else
    projected = isolated_refs_.first->coeff(); //so far, members of isolated_refs_ are the same

  // orthonormalize the "projected" coefficients
  Overlap S(sgeom_);
  Matrix S_invhalf = (*projected) % S * (*projected);
  S_invhalf.inverse_half();

  return make_shared<Matrix>(*projected * S_invhalf);
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
    Amatrix->copy_block(0, 0, dimerbasis, nclosed, scoeff->slice(0, nclosed)); // Total closed space
    Amatrix->copy_block(0, nclosed, dimerbasis, filled_activeB, scoeff->slice(nclosed+nactA, nclosed+nactA+filled_activeB)); // FilledActive B
    Amatrix->copy_block(0, nclosed+filled_activeB, dimerbasis, nactA, scoeff->slice(nclosed, nclosed+nactA)); // Active A
    auto Acoeff = make_shared<Coeff>(move(*Amatrix));

    // Set up variables for this fci
    const int ncore = nclosed + filled_activeB;
    const int norb  = nactA;

    embedded_refs_.first = make_shared<Reference>(sgeom_, Acoeff, ncore, norb, 0);
  }

  { // Move occupied orbitals of unit A to form core of unit B
    auto Bmatrix = make_shared<Matrix>(dimerbasis, nclosed + filled_activeA + nactB);
    Bmatrix->copy_block(0, 0, dimerbasis, nclosed, scoeff->slice(0, nclosed)); // Total closed space
    Bmatrix->copy_block(0, nclosed, dimerbasis, filled_activeA, scoeff->slice(nclosed, nclosed+filled_activeA)); // FilledActive A
    Bmatrix->copy_block(0, nclosed+filled_activeA, dimerbasis, nactB, scoeff->slice(nclosed+nactA, nclosed+nactA+nactB)); // Active B
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
    for (auto& s : *space)
      spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")});
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    for (auto& s : *spacea)
      spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")});
    for (auto& s : *spaceb)
      spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")});
  }
}


shared_ptr<Reference> Dimer::build_reference(const int site, const vector<bool> meanfield) const {
  const int nsites = meanfield.size();
  assert(nsites==2 && (site==0 || site==1));

  vector<shared_ptr<const MatView>> closed_orbitals = {make_shared<MatView>(sref_->coeff()->slice(0, sref_->nclosed()))};
  const MatView active_orbitals = sref_->coeff()->slice(sref_->nclosed() + (site==0 ? 0 : active_refs_.first->nact()),
                                                        sref_->nclosed() + active_refs_.first->nact() + (site==0 ? 0 : active_refs_.second->nact()));

  int current = sref_->nclosed();
  for (int i = 0; i < nsites; ++i) {
    const int cur_nact = (i==0 ? active_refs_.first->nact() : active_refs_.second->nact());
    const int cur_nocc = (i==0 ? isolated_refs_.first->nclosed() - active_refs_.first->nclosed() : isolated_refs_.first->nclosed() - active_refs_.second->nclosed());

    if (i!=site && meanfield[i])
      closed_orbitals.push_back(make_shared<const MatView>(sref_->coeff()->slice(current, current+cur_nocc)));
    current += cur_nact;
  }

  const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [] (const int a, shared_ptr<const MatView> m) { return a + m->mdim(); });
  const int nact = active_orbitals.mdim();

  auto out = make_shared<Matrix>(sref_->geom()->nbasis(), nclosed+nact);

  current = 0;
  closed_orbitals.push_back(make_shared<MatView>(active_orbitals));
  for (auto& orbitals : closed_orbitals) {
    copy_n(orbitals->data(), orbitals->mdim()*orbitals->ndim(), out->element_ptr(0, current));
    current += orbitals->mdim();
  }

  return make_shared<Reference>(sgeom_, make_shared<Coeff>(move(*out)), nclosed, nact, 0);
}


void Dimer::update_coeff(shared_ptr<const Matrix> new_coeff) {
  shared_ptr<Reference> temp_ref;
  temp_ref = make_shared<Reference>(*sref_);
  temp_ref->set_coeff(new_coeff);
  sref_ = make_shared<const Reference>(*temp_ref);
}
