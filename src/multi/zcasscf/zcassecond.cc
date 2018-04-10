//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcassecond.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/multi/zcasscf/zcassecond.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/mat1e/giao/reloverlap_london.h>
#include <src/util/math/quatmatrix.h>

using namespace std;
using namespace bagel;


ZCASSecond_base::ZCASSecond_base(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
 : ZCASSCF(idat, geom, ref) {
  thresh_microstep_ = idata_->get<double>("thresh_microstep", 1.0e-4);
  dfpcmo_ = idata_->get<bool>("dfpcmo", false);
}


ZCASSecond::ZCASSecond(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
 : ZCASSecond_base(idat, geom, ref) {
  init();
  cout << "   * Using the second-order algorithm" << endl << endl;
}


ZCASSecond_London::ZCASSecond_London(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
 : ZCASSecond_base(idat, geom, ref) {
  init();
  cout << "   * Using the second-order algorithm" << endl;
  cout << "   * A magnetic field is applied" << endl << endl;
}


void ZCASSecond::init_mat1e() {
  hcore_   = make_shared<RelHcore>(geom_);
  overlap_ = make_shared<RelOverlap>(geom_);
}


void ZCASSecond_London::init_mat1e() {
  hcore_   = make_shared<RelHcore_London>(geom_);
  overlap_ = make_shared<RelOverlap_London>(geom_);
}


void ZCASSecond::init_coeff() {
  auto relref = dynamic_pointer_cast<const RelReference>(ref_);

  // local function
  auto remove_lindep = [this](const int ndel) {
    assert(ndel > 0);
    nbasis_ -= 2*ndel;
    nneg_   -= 2*ndel;
    nvirt_  -= 2*ndel;
    nvirtnr_ -= ndel;
  };

  const bool hcore_guess = idata_->get<bool>("hcore_guess", false);
  shared_ptr<const ZCoeff_Striped> scoeff;
  if (hcore_guess) {
    auto s12 = overlap_->tildex(thresh_overlap_);
    if (s12->mdim() != s12->ndim())
      remove_lindep(geom_->nbasis() - s12->mdim()/4);

    auto hctmp = make_shared<ZMatrix>(*s12 % *hcore_ * *s12);
    VectorB eig(hctmp->ndim());
    hctmp->diagonalize(eig);
    scoeff = make_shared<const ZCoeff_Striped>(*s12 * *hctmp, nclosed_, nact_, nvirtnr_, nneg_, /*move_neg*/true);
  } else if (nr_coeff_ == nullptr) {

    // first set coefficient
    scoeff = relref->relcoeff_full();

    // If we have fewer MOs than AOs (linearly dependent basis set and/or coefficient projection)
    if (4 * geom_->nbasis() != scoeff->mdim()) {
      shared_ptr<const ZMatrix> tildex = overlap_->tildex(thresh_overlap_);
      if (tildex->mdim() != tildex->ndim())
        remove_lindep(geom_->nbasis() - tildex->mdim()/4);

      // Still not enough MOs (such as after coefficient projection)
      if (scoeff->mdim() != 2 * nbasis_) {
        ZMatrix c(scoeff->ndim(), tildex->mdim());
        c.copy_block(0, 0, scoeff->ndim(), scoeff->npos(), scoeff->slice(0, scoeff->npos()));
        c.copy_block(0, tildex->mdim() - scoeff->nneg(), scoeff->ndim(), scoeff->nneg(), scoeff->slice(scoeff->npos(), scoeff->mdim()));

        shared_ptr<const ZMatrix> trans = get<0>((*tildex % *overlap_ * *scoeff).svd());
        c.copy_block(0, scoeff->mdim(), scoeff->ndim(), tildex->mdim()-scoeff->mdim(), *tildex * trans->slice(scoeff->mdim(), tildex->mdim()));

        scoeff = make_shared<ZCoeff_Striped>(move(c), nclosed_, nact_, nvirtnr_, nneg_);
        scoeff = scoeff->init_kramers_coeff(geom_, overlap_, hcore_, geom_->nele() - charge_, gaunt_, breit_);
#ifndef NDEBUG
        ZMatrix unit(scoeff->mdim(), scoeff->mdim()); unit.unit();
        assert((*scoeff % *overlap_ * *scoeff - unit).rms() < 1.0e-10);
#endif
      }
    }

    scoeff = make_shared<const ZCoeff_Striped>(*scoeff, nclosed_, nact_, nvirtnr_, nneg_);
  } else {
    if (nr_coeff_->ndim() != nr_coeff_->mdim())
      remove_lindep(geom_->nbasis() - nr_coeff_->mdim());

    scoeff = nonrel_to_relcoeff(nr_coeff_)->striped_format();
  }

  // initialize coefficient to enforce kramers symmetry
  if (!relref->kramers() && (external_rdm_.empty() || external_rdm_ == "noref"))
    scoeff = scoeff->init_kramers_coeff(geom_, overlap_, hcore_, 2*ref_->nclosed() + ref_->nact(), gaunt_, breit_);

  // specify active orbitals and move into the active space
  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_);
  }
  coeff_ = scoeff->block_format();
}


void ZCASSecond_London::init_coeff() {
  auto relref = dynamic_pointer_cast<const RelReference>(ref_);

  // TODO this function is not as robust as it should be
  if (nr_coeff_)
    throw runtime_error("ZCASSecond_London requires a relativistic reference");

  // local function
  auto remove_lindep = [this](const int ndel) {
    assert(ndel > 0);
    nbasis_ -= 2*ndel;
    nneg_   -= 2*ndel;
    nvirt_  -= 2*ndel;
    nvirtnr_ -= ndel;
  };

  const bool hcore_guess = idata_->get<bool>("hcore_guess", false);
  shared_ptr<const ZCoeff_Striped> scoeff;
  if (hcore_guess) {
    auto s12 = overlap_->tildex(thresh_overlap_);
    if (s12->mdim() != s12->ndim())
      remove_lindep(geom_->nbasis() - s12->mdim()/4);
    auto hctmp = make_shared<ZMatrix>(*s12 % *hcore_ * *s12);
    VectorB eig(hctmp->ndim());
    hctmp->diagonalize(eig);
    scoeff = make_shared<const ZCoeff_Striped>(*s12 * *hctmp, nclosed_, nact_, nvirtnr_, nneg_, /*move_neg*/true);
  } else {
    scoeff = relref->relcoeff_full();
    scoeff = make_shared<const ZCoeff_Striped>(*scoeff, nclosed_, nact_, nvirtnr_, nneg_);
  }

  // specify active orbitals and move into the active space
  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_);
  }
  coeff_ = scoeff->block_format();
}


void ZCASSecond::impose_symmetry(shared_ptr<ZMatrix> o) const {
  assert(o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt_)*2 == o->ndim());

  auto kramers_adapt_block = [this,&o](const int jst, const int ist, const int joff, const int ioff) {
    for (int i = 0; i != ist; ++i)
      for (int j = 0; j != jst; ++j) {
        o->element(joff+j, ioff+i) = 0.5*(o->element(joff+j, ioff+i) + conj(o->element(joff+j+jst, ioff+i+ist)));
        o->element(joff+j+jst,ioff+i+ist) = conj(o->element(joff+j,ioff+i));

        o->element(joff+jst+j, ioff+i) = 0.5*(o->element(joff+jst+j, ioff+i) - conj(o->element(joff+j, ioff+ist+i)));
        o->element(joff+j, ioff+ist+i) = - conj(o->element(joff+jst+j, ioff+i));
      }
  };

  const array<int,3> stride{{nclosed_, nact_, nvirt_}};
  const array<int,3> offset{{0, 2*nclosed_, 2*nocc_}};
  for (int ii = 0; ii != 3; ++ii)
    for (int jj = 0; jj != 3; ++jj)
      kramers_adapt_block(stride[jj], stride[ii], offset[jj], offset[ii]);
}


void ZCASSecond::impose_symmetry(shared_ptr<ZRotFile> o) const {
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_; ++j) {
      o->ele_vc(j, i) = (o->ele_vc(j, i) + conj(o->ele_vc(j+nvirt_, i+nclosed_))) * 0.5;
      o->ele_vc(j+nvirt_, i+nclosed_) = conj(o->ele_vc(j, i));

      o->ele_vc(j+nvirt_, i) = (o->ele_vc(j+nvirt_, i) - conj(o->ele_vc(j, i+nclosed_))) * 0.5;
      o->ele_vc(j, i+nclosed_) = - conj(o->ele_vc(j+nvirt_, i));
    }
  }
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_; ++j) {
      o->ele_va(j, i) = (o->ele_va(j, i) + conj(o->ele_va(j+nvirt_, i+nact_))) * 0.5;
      o->ele_va(j+nvirt_, i+nact_) = conj(o->ele_va(j, i));

      o->ele_va(j+nvirt_, i) = (o->ele_va(j+nvirt_, i) - conj(o->ele_va(j, i+nact_))) * 0.5;
      o->ele_va(j, i+nact_) = - conj(o->ele_va(j+nvirt_, i));
    }
  }
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nclosed_; ++j) {
      o->ele_ca(j, i) = (o->ele_ca(j, i) + conj(o->ele_ca(j+nclosed_, i+nact_))) * 0.5;
      o->ele_ca(j+nclosed_, i+nact_) = conj(o->ele_ca(j, i));

      o->ele_ca(j+nclosed_, i) = (o->ele_ca(j+nclosed_, i) - conj(o->ele_ca(j, i+nact_))) * 0.5;
      o->ele_ca(j, i+nact_) = - conj(o->ele_ca(j+nclosed_, i));
    }
  }
}


