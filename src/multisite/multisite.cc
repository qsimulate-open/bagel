//
// BAGEL - Parallel electron correlation program.
// Filename: multisite.cc
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
#include <src/molecule/overlap.h>
#include <src/util/muffle.h>

using namespace std;
using namespace bagel;

MultiSite::MultiSite(shared_ptr<const PTree> input, vector<shared_ptr<const Reference>> refs) : input_(input), nsites_(refs.size()) {
  // build the super geometry
  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  vector<shared_ptr<const Geometry>> geovec;
  for (auto& r : refs) {
    geovec.push_back(mdata ? make_shared<Geometry>(*r->geom(), mdata) : r->geom());
    geoms_.push_back(r->geom());
  }

  auto envdata = input_->get_child_optional("environment");
  if (envdata) {
    Muffle hide_cout;
    geovec.push_back(make_shared<Geometry>(envdata));
  }

  auto sgeom = make_shared<Geometry>(geovec);

  // build combined references
  for (auto& r : refs)
    isolated_refs_.push_back(r->project_coeff(sgeom));
  const size_t norb = accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0ul, [] (size_t x, shared_ptr<const Reference> r) { return x + r->coeff()->mdim(); });
  Matrix coeff(sgeom->nbasis(), norb);

  // pull out all of the orbitals and put them in closed, active, virtual order
  int current = 0;
  for (auto& r : isolated_refs_) {
    copy_n(r->coeff()->element_ptr(0,0),                      r->nclosed()*sgeom->nbasis(), coeff.element_ptr(0, current));
    closed_bounds_.push_back({current, current+r->nclosed()});
    current += r->nclosed();
  }
  for (auto& r : isolated_refs_) {
    copy_n(r->coeff()->element_ptr(0,r->nclosed()+r->nact()), r->nvirt()*sgeom->nbasis(),   coeff.element_ptr(0, current));
    virt_bounds_.push_back({current, current+r->nvirt()});
    current += r->nvirt();
  }

  Overlap S(sgeom);
  Matrix transform(coeff % S * coeff);
  transform.inverse_half();
  auto scoeff = make_shared<Coeff>(coeff * transform);

  sref_ = make_shared<Reference>(sgeom, scoeff,
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nclosed(); }),
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nact(); }),
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nvirt(); }));
}

#if 0 // I'll do the SCF stuff first
shared_ptr<Reference> Dimer::build_reference(const int site, const vector<bool> meanfield) const {
  const int nsites = meanfield.size();
  assert(nsites==2 && (site==0 || site==1));

  vector<shared_ptr<const MatView>> closed_orbitals = {make_shared<MatView>(sref_->coeff()->slice(0, sref_->nclosed()))};
  const MatView active_orbitals = sref_->coeff()->slice(sref_->nclosed() + (site==0 ? 0 : active_refs_.first->nact()), sref_->nclosed() + (site==0 ? active_refs_.first->nact() : active_refs_.first->nact() + active_refs_.second->nact()));

  int current = sref_->nclosed();
  for (int i = 0; i < nsites; ++i) {
    const int cur_nact = (i==0 ? active_refs_.first->nact() : active_refs_.second->nact());
    const int cur_nocc = (i==0 ? isolated_refs_.first->nclosed() - active_refs_.first->nact() : isolated_refs_.first->nclosed() - active_refs_.second->nclosed());

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
#endif
