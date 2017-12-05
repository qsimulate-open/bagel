//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
// Copyright (C) 2014 Shane Parker
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

#include <src/asd/multisite/multisite.h>
#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;

///> this constructor takes localized HF orbitals with a clear definition of active subspaces
MultiSite::MultiSite(shared_ptr<const PTree> input, shared_ptr<const Reference> ref, const int nsites) : input_(input), hf_ref_(ref), nsites_(nsites) {
  cout << string(60, '=') << endl;
  cout << string(15, ' ') << "Construct MultiSite" << endl;
  cout << string(60, '=') << endl;

  charge_ = input_->get<int>("charge", 0);
  nspin_ = input_->get<int>("spin", 0);

  active_electrons_ = input_->get_vector<int>("active_electrons");
  active_sizes_ = input_->get_vector<int>("active_sizes");

  const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
  const int nclosed = (hf_ref_->geom()->nele() - charge_ - nactele) / 2;
  assert((hf_ref_->geom()->nele() - charge_ - nactele) % 2 == 0);

  const int nactive = accumulate(active_sizes_.begin(), active_sizes_.end(), 0);
  
  const int nvirt = hf_ref_->geom()->nbasis() - nclosed - nactive;

  sref_ = make_shared<Reference>(hf_ref_->geom(), nullptr, nclosed, nactive, nvirt);
}

void MultiSite::compute() {
  // construct Fock Matrix for canonicalization before setting active
  shared_ptr<const Matrix> density = hf_ref_->coeff()->form_density_rhf(sref_->nclosed());
  const MatView ccoeff = hf_ref_->coeff()->slice(0, hf_ref_->nclosed());
  auto fock = make_shared<const Fock<1>>(hf_ref_->geom(), hf_ref_->hcore(), density, ccoeff);

  // reorder coefficient for ASD-DMRG
  set_active();

  canonicalize(fock);

  run_fci();
}


void MultiSite::set_active() {
  auto active = input_->get_child("active");
  if (!active)
    throw logic_error("\"active\" keyword must be specified for multisite construction");
  if (active->size() != nsites_)
    throw logic_error("Must specify active spaces for all of the sites");

  vector<set<int>> active_orbitals;
  set<int> active_set;
  for (auto site : *active) {
    vector<int> active_orbs = site->get_vector<int>("");
    for_each(active_orbs.begin(), active_orbs.end(), [](int& x) { x--; });
    active_orbitals.emplace_back(active_orbs.begin(), active_orbs.end());
    active_set.insert(active_orbs.begin(), active_orbs.end());
  }
  
  const int nclosed = sref_->nclosed();
  const int nactive = sref_->nact();
  int closed_position = 0;
  int active_position = nclosed;
  int virt_position = nclosed + nactive;
  const int multisitebasis = hf_ref_->geom()->nbasis();

  auto hf_coeff = hf_ref_->coeff();
  auto out_coeff = hf_coeff->clone();
  assert(hf_coeff->ndim() == multisitebasis);

  for (auto site : active_orbitals) {
    for (int aorb : site)
      copy_n(hf_coeff->element_ptr(0, aorb), multisitebasis, out_coeff->element_ptr(0, active_position++));
  }
  for (int i = 0; i != hf_coeff->mdim(); ++i) {
    if (active_set.count(i) == 0)
      copy_n(hf_coeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, (closed_position < nclosed ? closed_position++ : virt_position++)));
  }
  assert(virt_position == hf_coeff->mdim());

  sref_ = make_shared<Reference>(hf_ref_->geom(), make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, sref_->nvirt());
}


void MultiSite::canonicalize(shared_ptr<const Matrix> fock) {
  //canonicalize active orbitals within subspace
  const int nclosed = sref_->nclosed();
  
  // rotate within each active subspace
  MatView active_mos = sref_->coeff()->slice(nclosed, nclosed+sref_->nact());
  Matrix fock_mo(active_mos % *fock * active_mos);
  VectorB eigs(active_mos.mdim());
  shared_ptr<Matrix> active_transformation = fock_mo.diagonalize_blocks(eigs, active_sizes_);
  Matrix transformed_mos(active_mos * *active_transformation);
  shared_ptr<Matrix> scoeff = sref_->coeff()->copy();
  scoeff->copy_block(0, nclosed, scoeff->ndim(), transformed_mos.mdim(), transformed_mos);

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


#include <src/util/muffle.h>
#include <src/ci/fci/knowles.h>
void MultiSite::run_fci() const {
  auto fci_info = input_->get_child_optional("fci");
  if(!fci_info) {
    cout << "No FCI info provided, skipping MultiSite debugging.." << endl;
    return;
  }

  auto fci = make_shared<KnowlesHandy>(fci_info, sref_->geom(), sref_);
  Muffle hide_cout("fci.log");
  fci->compute();
  hide_cout.unmute();
  cout << " * FCI energy : " << setprecision(12) << fci->energy(0) << endl;
}
