//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcoeff.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/wfn/zcoeff.h>
#include <src/util/math/quatmatrix.h>
#include <src/scf/dhf/dfock.h>
#include <cassert>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace bagel;

ZCoeff_base::ZCoeff_base(const int ndim, const bool loc, const int nclosed, const int nact, const int nvirt, const int nneg)
 : ZMatrix(ndim, 2*(nclosed+nact+nvirt)+nneg, loc), nbasis_(ndim/4), nclosed_(nclosed), nact_(nact), nvirt_nr_(nvirt), nneg_(nneg) {
  assert(ndim%4 == 0);
  assert(nneg%2 == 0);
  assert(npos() == nneg || nneg == 0 || nvirt_nr_ == 0);
}


void ZCoeff_base::print_info() const {
  cout << "    * nbasis   :  4 x " << nbasis_nr() << " = " << nbasis_rel() << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed() << endl;
  cout << "    * nact     : " << setw(6) << nact() << endl;
  cout << "    * nvirt_nr : " << setw(6) << nvirt_nr() << endl;
  cout << "    * nneg     : " << setw(6) << nneg() << endl;
}


ZCoeff_Striped::ZCoeff_Striped(const ZMatView& coeff, const int nclosed, const int nact, const int nvirt, const int nneg, const bool move_neg)
 : ZCoeff_base(coeff.ndim(), coeff.localized(), nclosed, nact, nvirt, nneg) {
  assert(ndim() == coeff.ndim());
  assert(mdim() == coeff.mdim());
  if (!move_neg) {
    // copy input matrix directly
    copy_block(0, 0, ndim(), mdim(), coeff);
  } else {
    // move positronic orbitals to end of virtual space
    copy_block(0, 0,      ndim(), npos(), coeff.element_ptr(0, nneg_));
    copy_block(0, npos(), ndim(), nneg_,  coeff.element_ptr(0, 0));
  }
}


ZCoeff_Block::ZCoeff_Block(const ZMatView& coeff, const int nclosed, const int nact, const int nvirt, const int nneg)
 : ZCoeff_base(coeff.ndim(), coeff.localized(), nclosed, nact, nvirt, nneg) {
  copy_block(0, 0, ndim(), mdim(), coeff);
}


ZCoeff_Kramers::ZCoeff_Kramers(const ZMatView& coeff, const int nclosed, const int nact, const int nvirt, const int nneg, const bool move_neg)
 : ZCoeff_base(coeff.ndim(), coeff.localized(), nclosed, nact, nvirt, nneg) {
  assert(nneg == npos());
  // copy input matrix directly
  copy_block(0, 0, ndim(), mdim(), coeff);
}


shared_ptr<ZCoeff_Block> ZCoeff_Striped::block_format(int nclosed, int nact, int nvirt_nr, int nneg) const {
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
  auto out = make_shared<ZCoeff_Block>(ndim(), localized(), nclosed, nact, nvirt_nr, nneg);
  // closed
  for (int j = 0; j != nclosed; ++j) {
    out->copy_block(0,           j, n, 1, slice(j*2  , j*2+1));
    out->copy_block(0, nclosed + j, n, 1, slice(j*2+1, j*2+2));
  }
  int offset = nclosed*2;
  // active
  for (int j = 0; j != nact; ++j) {
    out->copy_block(0, offset + j,        n, 1, slice(offset +j*2,   offset + j*2+1));
    out->copy_block(0, offset + nact + j, n, 1, slice(offset +j*2+1, offset + j*2+2));
  }
  offset = (nclosed+nact)*2;
  // virtual (including positrons)
  for (int j = 0; j != nvirt_nr + nneg/2; ++j) {
    out->copy_block(0, offset + j,                     n, 1, slice(offset + j*2,   offset + j*2+1));
    out->copy_block(0, offset + nvirt_nr+nneg/2 + j,   n, 1, slice(offset + j*2+1, offset + j*2+2));
  }
  return out;
}


shared_ptr<ZCoeff_Striped> ZCoeff_Block::striped_format() const {
  assert(nneg_ % 2 == 0);
  int n = ndim();
  int offset = nclosed_;
  auto out = make_shared<ZCoeff_Striped>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);
  // closed
  for (int j = 0; j != nclosed_; ++j) {
    out->copy_block(0, j*2,   n, 1, slice(j, j+1));
    out->copy_block(0, j*2+1, n, 1, slice(offset + j, offset + j+1));
  }
  offset = nclosed_*2;
  // active
  for (int j = 0; j != nact_; ++j) {
    out->copy_block(0, offset + j*2,   n, 1, slice(offset + j,         offset + j+1));
    out->copy_block(0, offset + j*2+1, n, 1, slice(offset + nact_ + j, offset + nact_ + j+1));
  }
  offset = (nclosed_+nact_)*2;
  // vituals (including positrons)
  for (int j = 0; j != nvirt_rel(); ++j) {
    out->copy_block(0, offset + j*2,   n, 1, slice(offset + j,          offset + j+1));
    out->copy_block(0, offset + j*2+1, n, 1, slice(offset + nvirt_rel() + j, offset + nvirt_rel() + j+1));
  }
  return out;
}


shared_ptr<ZCoeff_Striped> ZCoeff_Kramers::striped_format() const {
  assert(nneg_ % 2 == 0);
  int n = ndim();
  int offset = nclosed_ + nact_ + nvirt_nr_ + nneg_/2;
  auto out = make_shared<ZCoeff_Striped>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);

  for (int j = 0; j != offset; ++j) {
    out->copy_block(0, j*2,   n, 1, slice(j, j+1));
    out->copy_block(0, j*2+1, n, 1, slice(offset + j, offset + j+1));
  }
  return out;
}


shared_ptr<ZCoeff_Block> ZCoeff_Kramers::block_format() const {
  int i = 0;
  array<shared_ptr<const ZMatrix>,2> tmp = {{ slice_copy(0, mdim()/2), slice_copy(mdim()/2, mdim()) }};
  auto out = make_shared<ZCoeff_Block>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);

  out->copy_block(0, i, ndim(), nclosed_, tmp[0]->slice(0,nclosed_)); i += nclosed_;
  out->copy_block(0, i, ndim(), nclosed_, tmp[1]->slice(0,nclosed_)); i += nclosed_;
  out->copy_block(0, i, ndim(), nact_, tmp[0]->slice(nclosed_, nclosed_+nact_)); i += nact_;
  out->copy_block(0, i, ndim(), nact_, tmp[1]->slice(nclosed_, nclosed_+nact_)); i += nact_;
  out->copy_block(0, i, ndim(), nvirt_rel(), tmp[0]->slice(nclosed_+nact_, tmp[0]->mdim())); i += nvirt_rel();
  out->copy_block(0, i, ndim(), nvirt_rel(), tmp[1]->slice(nclosed_+nact_, tmp[1]->mdim()));
  return out;
}


shared_ptr<ZCoeff_Kramers> ZCoeff_Kramers::swap_central() const {
  auto out = make_shared<ZCoeff_Kramers>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);
  const int m = nclosed_ + nact_ + nvirt_nr_;

  out->copy_block(0, 0*m, ndim(), m, slice(0*m, 1*m));
  out->copy_block(0, 1*m, ndim(), m, slice(2*m, 3*m));
  out->copy_block(0, 2*m, ndim(), m, slice(1*m, 2*m));
  out->copy_block(0, 3*m, ndim(), m, slice(3*m, 4*m));
  return out;
}


shared_ptr<ZCoeff_Kramers> ZCoeff_Kramers::move_positronic() const {
  auto out = make_shared<ZCoeff_Kramers>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, nneg_);
  const int m = nclosed_ + nact_ + nvirt_nr_;

  out->copy_block(0, 0*m, ndim(), m, slice(1*m, 2*m));
  out->copy_block(0, 1*m, ndim(), m, slice(0*m, 1*m));
  out->copy_block(0, 2*m, ndim(), m, slice(3*m, 4*m));
  out->copy_block(0, 3*m, ndim(), m, slice(2*m, 3*m));
  return out;
}


shared_ptr<Kramers<1,ZMatrix>> ZCoeff_Striped::kramers_active() const {
  ZCoeff_Striped active_only(slice(2*nclosed_, 2*(nclosed_+nact_)), 0, nact_, 0, 0);
  return active_only.block_format()->kramers_active();
}


shared_ptr<Kramers<1,ZMatrix>> ZCoeff_Block::kramers_active() const {
  auto out = make_shared<Kramers<1,ZMatrix>>();
  out->emplace(0, slice_copy(2*nclosed_,         2*nclosed_ + nact_));
  out->emplace(1, slice_copy(2*nclosed_ + nact_, 2*(nclosed_+nact_)));
  return out;
}


shared_ptr<ZCoeff_Block> ZCoeff_Block::electronic_part() const {
  auto out = make_shared<ZCoeff_Block>(ndim(), localized(), nclosed_, nact_, nvirt_nr_, 0);
  out->copy_block(0, 0,                  ndim(), 2*nocc(),  slice(0,                    2*nocc()));
  out->copy_block(0, 2*nocc(),           ndim(), nvirt_nr_, slice(2*nocc(),             2*nocc()+nvirt_nr_));
  out->copy_block(0, 2*nocc()+nvirt_nr_, ndim(), nvirt_nr_, slice(2*nocc()+nvirt_rel(), 2*nocc()+nvirt_rel()+nvirt_nr_));
  return out;
}


shared_ptr<ZCoeff_Block> ZCoeff_Block::closed_part() const {
  auto out = make_shared<ZCoeff_Block>(slice(0, 2*nclosed_), nclosed_, 0, 0, 0);
  return out;
}


shared_ptr<ZCoeff_Block> ZCoeff_Block::active_part() const {
  auto out = make_shared<ZCoeff_Block>(slice(2*nclosed_, 2*nocc()), 0, nact_, 0, 0);
  return out;
}


// Kramers-adapted coefficient via quaternion diagonalization, assuming guess orbitals from Dirac--Hartree--Fock
shared_ptr<const ZCoeff_Striped> ZCoeff_Striped::init_kramers_coeff(shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> overlap,
                                                                    shared_ptr<const ZMatrix> hcore, const int nele, const bool gaunt, const bool breit) const {

  // quaternion diagonalization has a bug for 2x2 case since there are no super-offdiagonals in a 2x2 and tridiagonalization is probably not possible
  assert(nact_ != 1);
  assert(nbasis_ == geom->nbasis());

  shared_ptr<const ZMatrix> orthog = overlap->tildex(1.0e-9);
  shared_ptr<const ZMatrix> s12 = ZCoeff_Kramers(*orthog, nclosed_, nact_, nvirt_nr(), nneg_).swap_central();
  shared_ptr<const ZMatrix> focktmp = make_shared<DFock>(geom, hcore, slice_copy(0, nele), gaunt, breit, /*store_half*/false, /*robust*/breit);
  VectorB eig(s12->mdim());

  // quaternion diagonalize a fock matrix in MO basis
  shared_ptr<ZMatrix> fock_tilde;
  {
    fock_tilde = make_shared<QuatMatrix>(*s12 % (*focktmp) * *s12);
#ifndef NDEBUG
    auto quatfock = static_pointer_cast<const QuatMatrix>(fock_tilde);
    const double tsymm_err = quatfock->check_t_symmetry();
    if (tsymm_err > 1.0e-8) {
      cout << "   ** Caution:  poor Kramers symmetry in Fock matrix used to symmetrize orbitals:  error = " << scientific << setprecision(4) << tsymm_err << endl;
      cout << "   ** This may result in poor quality orbitals for ZFCI or the ZCASSCF initial guess." << endl;
    }
#endif
  }

  fock_tilde->diagonalize(eig);

  // move positronic orbitals to the end and transform to striped format
  auto out = make_shared<const ZCoeff_Kramers>(*s12 * *fock_tilde, nclosed_, nact_, nvirt_nr_, nneg_);
  return out->move_positronic()->striped_format();
}


shared_ptr<const ZCoeff_Striped> ZCoeff_Striped::set_active(set<int> active_indices, const int nele) const {
  const int nmobasis = npos()/2;

  cout << " " << endl;
  cout << "    ==== Active orbitals : ===== " << endl;
  for (auto& i : active_indices) cout << "         Orbital " << i+1 << endl;
  cout << "    ============================ " << endl << endl;

  if (active_indices.size() != nact_)
    throw logic_error("ZCoeff_Striped::set_active - Number of active indices does not match number of active orbitals.  (" + to_string(nact_) + " expected)");
  if (any_of(active_indices.begin(), active_indices.end(), [nmobasis](int i){ return (i < 0 || i >= nmobasis); }) )
    throw runtime_error("ZCoeff_Striped::set_active - Invalid MO index provided.  (Should be from 1 to " + to_string(nmobasis) + ")");

  auto out = make_shared<ZCoeff_Striped>(ndim(), localized(), nclosed_, nact_, mdim()/4-nclosed_-nact_, nneg());

  int iclosed = 0;
  int iactive = nclosed_;
  int ivirt   = nclosed_ + nact_;

  auto cp = [&out, this] (const int i, int& pos) {
    copy_n(element_ptr(0,i*2), nbasis_rel(), out->element_ptr(0, pos*2));
    copy_n(element_ptr(0,i*2+1), nbasis_rel(), out->element_ptr(0, pos*2+1));
    ++pos;
  };

  int closed_count = 0;
  for (int i = 0; i < nmobasis; ++i) {
    if (active_indices.find(i) != active_indices.end()) {
      cp(i, iactive);
    } else if (closed_count < nclosed_) {
      cp(i, iclosed);
      closed_count++;
    } else {
      cp(i, ivirt);
    }
  }

  if (closed_count != nclosed_)
    throw runtime_error("Invalid combination of closed and active orbitals.");

  // copy positrons
  out->copy_block(0, npos(), nbasis_rel(), nneg(), slice(npos(), mdim()));
  return out;
}


void ZCoeff_base::rearrange_eig(VectorB& eig, shared_ptr<ZMatrix> coeff, const bool includes_neg) {
  const int n = coeff->ndim()/2;
  assert(2*n == coeff->ndim());  // could be triggered if Kramers + and - sets had different sizes or linear dependencies

  // check that pos & neg energy eigenvalues are properly separated
  assert(!includes_neg || *std::min_element(eig.begin()+n, eig.begin()+2*n) - *std::max_element(eig.begin(), eig.begin()+n) > c__*c__);

  // need to reorder things so negative energy states don't all come at the beginning
  VectorB tempv(2*n);
  shared_ptr<ZMatrix> tempm = coeff->clone();
  for (int i = 0; i != n; ++i) {
    tempv[  i] = eig[2*i];
    tempv[n+i] = eig[2*i+1];
    tempm->copy_block(0,   i, 2*n, 1, coeff->slice(2*i,   2*i+1));
    tempm->copy_block(0, n+i, 2*n, 1, coeff->slice(2*i+1, 2*i+2));
  }
  eig = tempv;
  *coeff = *tempm;
}


BOOST_CLASS_EXPORT_IMPLEMENT(ZCoeff_base)
BOOST_CLASS_EXPORT_IMPLEMENT(ZCoeff_Striped)
BOOST_CLASS_EXPORT_IMPLEMENT(ZCoeff_Block)
BOOST_CLASS_EXPORT_IMPLEMENT(ZCoeff_Kramers)
