//
// BAGEL - Parallel electron correlation program.
// Filename: relcoeff.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/wfn/relcoeff.h>
#include <src/util/math/quatmatrix.h>
#include <src/scf/dhf/dfock.h>
#include <cassert>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace bagel;

RelCoeff::RelCoeff(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
 : ZMatrix(_ndim, 2*(_nclosed+_nact+_nvirt)+_nneg, _loc), nbasis_(ndim()/4), nclosed_(_nclosed), nact_(_nact), nvirt_nr_(_nvirt), nneg_(_nneg) {
  assert(ndim()%4 == 0);
  assert(nneg()%2 == 0);
  assert(npos() == nneg() || nneg() == 0);
}


void RelCoeff::print_info() const {
  cout << "    * nbasis   :  4 x " << nbasis_nr() << " = " << nbasis_rel() << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed() << endl;
  cout << "    * nact     : " << setw(6) << nact() << endl;
  cout << "    * nvirt_nr : " << setw(6) << nvirt_nr() << endl;
  cout << "    * nneg     : " << setw(6) << nneg() << endl;
}


RelCoeff_Striped::RelCoeff_Striped(const ZMatView& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg, const bool move_neg)
 : RelCoeff(_coeff.ndim(), _coeff.localized(), _nclosed, _nact, _nvirt, _nneg) {
  if (!move_neg) {
    // copy input matrix directly
    copy_block(0, 0, ndim(), mdim(), _coeff);
  } else {
    // move positronic orbitals to end of virtual space
    copy_block(0, 0,      ndim(), npos(), _coeff.element_ptr(0, nneg_));
    copy_block(0, npos(), ndim(), nneg_,  _coeff.element_ptr(0, 0));
  }
}


RelCoeff_Block::RelCoeff_Block(const ZMatView& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
 : RelCoeff(_coeff.ndim(), _coeff.localized(), _nclosed, _nact, _nvirt, _nneg) {
  copy_block(0, 0, ndim(), mdim(), _coeff);
}


RelCoeff_Kramers::RelCoeff_Kramers(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
 : RelCoeff(_coeff.ndim(), _coeff.localized(), _nclosed, _nact, _nvirt, _nneg) {

  const int n = nbasis_nr();
  const int m = nclosed_ + nact_ + nvirt_nr_;
  assert(nneg_ == 0 || nneg_ == 2*m);

  auto copy_swapped_blocks = [this, &_coeff, n, m](const int ao1, const int mo1, const int ao2, const int mo2) {
    copy_block(ao1*n, mo1*m, n, m, _coeff.get_submatrix(ao2*n, mo2*m, n, m));
    copy_block(ao2*n, mo2*m, n, m, _coeff.get_submatrix(ao1*n, mo1*m, n, m));
  };

  //    rows: {L+, S+, L-, S-} -> {L+, L-, S+, S-}
  // columns: {S+, L+, S-, L-} -> {L+, S+, L-, S-}
  copy_swapped_blocks(0, 0, 0, 1);
  copy_swapped_blocks(0, 2, 0, 3);
  copy_swapped_blocks(1, 0, 2, 1);
  copy_swapped_blocks(1, 1, 2, 0);
  copy_swapped_blocks(1, 2, 2, 3);
  copy_swapped_blocks(1, 3, 2, 2);
  copy_swapped_blocks(3, 0, 3, 1);
  copy_swapped_blocks(3, 2, 3, 3);
}


std::shared_ptr<RelCoeff_Block> RelCoeff_Striped::block_format(int nclosed, int nact, int nvirt_nr, int nneg) const {
  if (nneg == -1) {
    assert(nclosed == -1 && nact == -1 && nvirt_nr == -1);
    nclosed = nclosed_;
    nact = nact_;
    nvirt_nr = nvirt_nr_;
    nneg = nneg_;
  }
  assert(nneg % 2 == 0);
  assert(2*(nclosed+nact+nvirt_nr) == nneg || nneg == 0);
  assert(2*(nclosed+nact+nvirt_nr) + nneg == mdim());

  const int n = ndim();
  auto out = make_shared<RelCoeff_Block>(ndim(), localized(), nclosed, nact, nvirt_nr, nneg);
  // closed
  for (int j=0; j!=nclosed; ++j) {
    out->copy_block(0,           j, n, 1, slice(j*2  , j*2+1));
    out->copy_block(0, nclosed + j, n, 1, slice(j*2+1, j*2+2));
  }
  int offset = nclosed*2;
  // active
  for (int j=0; j!=nact; ++j) {
    out->copy_block(0, offset + j,        n, 1, slice(offset +j*2,   offset + j*2+1));
    out->copy_block(0, offset + nact + j, n, 1, slice(offset +j*2+1, offset + j*2+2));
  }
  offset = (nclosed+nact)*2;
  // virtual (including positrons)
  for (int j=0; j!=nvirt_nr+nneg/2; ++j) {
    out->copy_block(0, offset + j,                     n, 1, slice(offset + j*2,   offset + j*2+1));
    out->copy_block(0, offset + nvirt_nr+nneg/2 + j,   n, 1, slice(offset + j*2+1, offset + j*2+2));
  }
  return out;
}


std::shared_ptr<RelCoeff_Striped> RelCoeff_Block::striped_format() const {
  assert(nneg_ % 2 == 0);
  int n = ndim();
  int offset = nclosed_;
  auto out = make_shared<RelCoeff_Striped>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);
  // closed
  for (int j=0; j!=nclosed_; ++j) {
    out->copy_block(0, j*2,   n, 1, slice(j, j+1));
    out->copy_block(0, j*2+1, n, 1, slice(offset + j, offset + j+1));
  }
  offset = nclosed_*2;
  // active
  for (int j=0; j!=nact_; ++j) {
    out->copy_block(0, offset + j*2,   n, 1, slice(offset + j,         offset + j+1));
    out->copy_block(0, offset + j*2+1, n, 1, slice(offset + nact_ + j, offset + nact_ + j+1));
  }
  offset = (nclosed_+nact_)*2;
  // vituals (including positrons)
  for (int j=0; j!=nvirt_rel(); ++j) {
    out->copy_block(0, offset + j*2,   n, 1, slice(offset + j,          offset + j+1));
    out->copy_block(0, offset + j*2+1, n, 1, slice(offset + nvirt_rel() + j, offset + nvirt_rel() + j+1));
  }
  return out;
}


std::shared_ptr<RelCoeff_Block> RelCoeff_Kramers::block_format() const {

  int i = 0;
  array<shared_ptr<const ZMatrix>,2> tmp = {{ slice_copy(0, mdim()/2), slice_copy(mdim()/2, mdim()) }};
  auto out = make_shared<RelCoeff_Block>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);

  out->copy_block(0, i, ndim(), nclosed_, tmp[0]->slice(0,nclosed_)); i += nclosed_;
  out->copy_block(0, i, ndim(), nclosed_, tmp[1]->slice(0,nclosed_)); i += nclosed_;
  out->copy_block(0, i, ndim(), nact_, tmp[0]->slice(nclosed_, nclosed_+nact_)); i += nact_;
  out->copy_block(0, i, ndim(), nact_, tmp[1]->slice(nclosed_, nclosed_+nact_)); i += nact_;
  out->copy_block(0, i, ndim(), nvirt_rel(), tmp[0]->slice(nclosed_+nact_, tmp[0]->mdim())); i += nvirt_rel();
  out->copy_block(0, i, ndim(), nvirt_rel(), tmp[1]->slice(nclosed_+nact_, tmp[1]->mdim()));

  return out;
}


std::shared_ptr<RelCoeff_Striped> RelCoeff_Kramers::striped_format() const {
  assert(nneg_ % 2 == 0);
  int n = ndim();
  int offset = nclosed_ + nact_ + nvirt_nr_ + nneg_/2;
  auto out = make_shared<RelCoeff_Striped>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);

  for (int j=0; j!=offset; ++j) {
    out->copy_block(0, j*2,   n, 1, slice(j, j+1));
    out->copy_block(0, j*2+1, n, 1, slice(offset + j, offset + j+1));
  }
  return out;
}


std::shared_ptr<Kramers<2,ZMatrix>> RelCoeff_Striped::kramers_active() const {
  RelCoeff_Striped active_only(slice(2*nclosed_, 2*(nclosed_+nact_)), 0, nact_, 0, 0);
  return active_only.block_format()->kramers_active();
}


std::shared_ptr<Kramers<2,ZMatrix>> RelCoeff_Block::kramers_active() const {
  auto out = make_shared<Kramers<2,ZMatrix>>();
  out->emplace(0, slice_copy(2*nclosed_,         2*nclosed_ + nact_));
  out->emplace(1, slice_copy(2*nclosed_ + nact_, 2*(nclosed_+nact_)));
  return out;
}


// Kramers-adapted coefficient via quaternion diagonalization, assuming guess orbitals from Dirac--Hartree--Fock
shared_ptr<const RelCoeff_Striped> RelCoeff_Striped::init_kramers_coeff_dirac(shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> overlap,
                          shared_ptr<const ZMatrix> hcore, const int nele, const bool tsymm, const bool gaunt, const bool breit) const {

  // quaternion diagonalization has a bug for 2x2 case since there are no super-offdiagonals in a 2x2 and tridiagonalization is probably not possible
  assert((nact_ > 1 || !tsymm));
  assert(nbasis_ == geom->nbasis());

  // transformation from the standard format to the quaternion format
  auto quaternion = [](shared_ptr<ZMatrix> o) {
    shared_ptr<ZMatrix> scratch = o->clone();
    const int n = o->ndim()/4;
    const int m = o->mdim()/4;
    map<int, int> trans {{0,0}, {1,2}, {2,1}, {3,3}};
    for (auto& i : trans)
      for (auto& j : trans)
        scratch->copy_block(i.first*n, j.first*m, n, m, o->get_submatrix(i.second*n, j.second*m, n, m));
    *o = *scratch;
  };

  shared_ptr<ZMatrix> coefftmp;
  shared_ptr<ZMatrix> focktmp;
  shared_ptr<ZMatrix> ctmp;
  // quaternion diagonalize a fock matrix
  focktmp = make_shared<DFock>(geom, hcore, slice_copy(0, nele), gaunt, breit, /*store_half*/false, /*robust*/false);
  quaternion(focktmp);

  shared_ptr<ZMatrix> s12 = overlap->tildex(1.0e-9);
  quaternion(s12);

  shared_ptr<ZMatrix> fock_tilde;
  if (tsymm) {
    fock_tilde = make_shared<QuatMatrix>(*s12 % (*focktmp) * *s12);
#ifndef NDEBUG
    auto quatfock = static_pointer_cast<const QuatMatrix>(fock_tilde);
    const double tsymm_err = quatfock->check_t_symmetry();
    if (tsymm_err > 1.0e-8)
      cout << "   ** Caution:  poor Kramers symmetry in fock_tilde (ZCASSCF initialization) - error = " << scientific << setprecision(4) << tsymm_err << endl;
#endif
  } else {
    fock_tilde = make_shared<ZMatrix>(*s12 % (*focktmp) * *s12);
  }

  // quaternion diagonalization
  {
    VectorB eig(fock_tilde->ndim());
    fock_tilde->diagonalize(eig);

    if (!tsymm)
      rearrange_eig(eig, fock_tilde);
  }

  // re-order to kramers format and move negative energy states to virtual space
  ctmp = make_shared<ZMatrix>(*s12 * *fock_tilde);

  auto ktmp = make_shared<const RelCoeff_Kramers>(*ctmp, nclosed_, nact_, nvirt_nr_, nneg_);
  auto out = ktmp->striped_format();
  return out;
}


shared_ptr<const RelCoeff_Striped> RelCoeff_Striped::set_active(set<int> active_indices, const int nele, const bool paired) const {
  // assumes coefficient is in striped format
  const int nmobasis = paired ? npos()/2 : npos();

  cout << " " << endl;
  if (!paired) cout << "    * Selecting individual spin-orbitals for the active space." << endl;
  cout << "    ==== Active orbitals : ===== " << endl;
  for (auto& i : active_indices) cout << "         Orbital " << i+1 << endl;
  cout << "    ============================ " << endl << endl;

  if (active_indices.size() != (paired ? nact_ : 2*nact_))
    throw logic_error("RelCoeff_Striped::set_active - Number of active indices does not match number of active orbitals.  (" + to_string(paired ? nact_ : 2*nact_) + " expected)");
  if (any_of(active_indices.begin(), active_indices.end(), [nmobasis](int i){ return (i < 0 || i >= nmobasis); }) )
    throw runtime_error("RelCoeff_Striped::set_active - Invalid MO index provided.  (Should be from 1 to " + to_string(nmobasis) + ")");

  shared_ptr<ZMatrix> tmp_coeff = clone();

  int iclosed = 0;
  int iactive = nclosed_;
  int ivirt   = nclosed_ + nact_;
  int nclosed_start = nele / 2;

  if (!paired) {
    iactive *= 2;
    ivirt *= 2;
    nclosed_start = nele;
  }

  auto cp   = [&tmp_coeff, this, &paired] (const int i, int& pos) {
    if (paired) {
      copy_n(element_ptr(0,i*2), nbasis_rel(), tmp_coeff->element_ptr(0, pos*2));
      copy_n(element_ptr(0,i*2+1), nbasis_rel(), tmp_coeff->element_ptr(0, pos*2+1));
    } else {
      copy_n(element_ptr(0,i), nbasis_rel(), tmp_coeff->element_ptr(0, pos));
    }
    ++pos;
  };

  int closed_count = 0;
  for (int i = 0; i < nmobasis; ++i) {
    if (active_indices.find(i) != active_indices.end()) {
      cp(i, iactive);
    } else if (i < nclosed_start) {
      cp(i, iclosed);
      closed_count++;
    } else {
      cp(i, ivirt);
    }
  }

  if (closed_count != (paired ? nclosed_ : 2*nclosed_))
    throw runtime_error("Invalid combination of closed and active orbitals.");

  // copy positrons
  tmp_coeff->copy_block(0, npos(), nbasis_rel(), nneg(), slice(npos(), mdim()));

  return make_shared<const RelCoeff_Striped>(*tmp_coeff, nclosed_, nact_, mdim()/4-nclosed_-nact_, nneg());
}


void RelCoeff::rearrange_eig(VectorB& eig, shared_ptr<ZMatrix> coeff, const bool includes_neg) {
  const int n = coeff->ndim()/2;
  assert(2*n == coeff->ndim());  // could be triggered if Kramers + and - sets had different sizes or linear dependencies

  // check that pos & neg energy eigenvalues are properly separated
  assert(!includes_neg || *std::min_element(eig.begin()+n, eig.begin()+2*n) - *std::max_element(eig.begin(), eig.begin()+n) > c__*c__);

  // need to reorder things so negative energy states don't all come at the beginning
  // TODO there should be a more efficient way...
  VectorB tempv(2*n);
  shared_ptr<ZMatrix> tempm = coeff->clone();
  for (int i = 0; i != n; ++i) {
    tempv[  i] = eig[2*i];
    tempv[n+i] = eig[2*i+1];

    for (int j=0; j!=2*n; ++j) {
      tempm->element(j,   i) = coeff->element(j, 2*i);
      tempm->element(j, n+i) = coeff->element(j, 2*i+1);
    }
  }
  eig = tempv;
  *coeff = *tempm;
}


// function to compute the modified virtual orbitals, either by diagonalization of a Fock matrix or of the one-electron Hamiltonian
// Procedures described in Jensen et al; JCP 87, 451 (1987) (hcore) and Bauschlicher; JCP 72 880 (1980) (Fock)
shared_ptr<const RelCoeff_Striped> RelCoeff_Striped::generate_mvo(shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> overlap,
        shared_ptr<const ZMatrix> hcore, const int ncore, const int nocc_mvo, const bool hcore_mvo, const bool tsymm, const bool gaunt, const bool breit) const {
  cout << " " << endl;
  if (!hcore_mvo) {
    cout << "   * Generating Modified Virtual Orbitals from a Fock matrix of " << ncore << " electrons " << endl << endl;
  } else {
    cout << "   * Generating Modified Virtual Orbitals from the 1 electron Hamiltonian of " << ncore << " electrons " << endl << endl;
  }
  assert(2*nocc_mvo >= ncore);
  const int hfvirt = nvirt_nr_ + nclosed_ + nact_ - nocc_mvo;
  assert(2*(nocc_mvo + hfvirt) + nneg_ == mdim());

  shared_ptr<const ZMatrix> mvofock = !hcore_mvo ? make_shared<const DFock>(geom, hcore, slice_copy(0, ncore), gaunt, breit, /*store half*/false, /*robust*/breit) : hcore;

  // take virtual part out and make block format
  shared_ptr<RelCoeff_Block> vcoeff = RelCoeff_Striped(slice(2*nocc_mvo, 2*(nocc_mvo + hfvirt)), 0, 0, hfvirt, 0).block_format();

  shared_ptr<ZMatrix> mofock;
  if (tsymm) {
    mofock = make_shared<QuatMatrix>(*vcoeff % *mvofock * *vcoeff);
#ifndef NDEBUG
    auto quatfock = static_pointer_cast<const QuatMatrix>(mofock);
    assert(quatfock->is_t_symmetric());
#endif
  } else {
    mofock = make_shared<ZMatrix>(*vcoeff % *mvofock * *vcoeff);
  }

  VectorB eig(mofock->ndim());
  mofock->diagonalize(eig);

  if (!tsymm)
    rearrange_eig(eig, mofock);

  // update orbitals and back transform
  shared_ptr<RelCoeff_Striped> scoeff = RelCoeff_Block(*vcoeff * *mofock, 0, 0, hfvirt, 0).striped_format();

  // copy in modified virtuals
  shared_ptr<ZMatrix> ecoeff = this->copy();
  ecoeff->copy_block(0, 2*nocc_mvo, ecoeff->ndim(), 2*hfvirt, scoeff->data());

  {
    auto unit = ecoeff->clone(); unit->unit();
    double orthonorm = ((*ecoeff % *overlap * *ecoeff) - *unit).rms();
    if (orthonorm > 1.0e-12) throw logic_error("MVO Coefficient not sufficiently orthonormal");
  }

  return make_shared<const RelCoeff_Striped>(*ecoeff, nclosed_, nact_, nvirt_nr_, nneg_);
}


BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Striped)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Block)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Kramers)
