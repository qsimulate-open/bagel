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


RelCoeff_Striped::RelCoeff_Striped(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg, const bool move_neg)
 : RelCoeff(_coeff.ndim(), _coeff.localized(), _nclosed, _nact, _nvirt, _nneg) {
  if (!move_neg) {
    // copy input matrix directly
    copy_block(0, 0, ndim(), mdim(), _coeff);
  } else {
    // move positronic orbitals to end of virtual space
    copy_block(0, 0,      ndim(), npos(), _coeff.slice(nneg_, nneg_+npos()));
    copy_block(0, npos(), ndim(), nneg_,  _coeff.slice(0,     nneg_));
  }
}


RelCoeff_Block::RelCoeff_Block(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
 : RelCoeff(_coeff.ndim(), _coeff.localized(), _nclosed, _nact, _nvirt, _nneg) {
  copy_block(0, 0, ndim(), mdim(), _coeff);
}


RelCoeff_Kramers::RelCoeff_Kramers(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
 : RelCoeff(_coeff.ndim(), _coeff.localized(), _nclosed, _nact, _nvirt, _nneg) {

  // TODO Should be a bit cheaper if we do all these moves together and get rid of scratch.
  // rows: {L+, S+, L-, S-} -> {L+, L-, S+, S-}
  copy_block(        0, 0, nbasis_, mdim(), _coeff.get_submatrix(        0, 0, nbasis_, mdim()));
  copy_block(  nbasis_, 0, nbasis_, mdim(), _coeff.get_submatrix(2*nbasis_, 0, nbasis_, mdim()));
  copy_block(2*nbasis_, 0, nbasis_, mdim(), _coeff.get_submatrix(  nbasis_, 0, nbasis_, mdim()));
  copy_block(3*nbasis_, 0, nbasis_, mdim(), _coeff.get_submatrix(3*nbasis_, 0, nbasis_, mdim()));

  // columns: {S+, L+, S-, L-} -> {L+, S+, L-, S-}
  auto move_one = [this](const int offset, const int block1, const int block2) {
    shared_ptr<ZMatrix> scratch = make_shared<ZMatrix>(ndim(), block1+block2);
    scratch->copy_block(0,      0, ndim(), block2, slice(offset+block1, offset+block1+block2));
    scratch->copy_block(0, block2, ndim(), block1, slice(offset,        offset+block1));
    copy_block(0, offset, ndim(), block1+block2, scratch);
  };
  move_one(       0, nneg()/2, npos()/2);
  move_one(mdim()/2, nneg()/2, npos()/2);
}


// Transforms a coefficient matrix from striped format to block format : assumes ordering is (c,a,v,positrons)
std::shared_ptr<RelCoeff_Block> RelCoeff_Striped::block_format() const {
  assert(nneg_ % 2 == 0);
  shared_ptr<ZMatrix> ctmp2 = clone();
  int n = ndim();
  // closed
  for (int j=0; j!=nclosed_; ++j) {
    ctmp2->copy_block(0,            j, n, 1, slice(j*2  , j*2+1));
    ctmp2->copy_block(0, nclosed_ + j, n, 1, slice(j*2+1, j*2+2));
  }
  int offset = nclosed_*2;
  // active
  for (int j=0; j!=nact_; ++j) {
    ctmp2->copy_block(0, offset + j,         n, 1, slice(offset +j*2,   offset + j*2+1));
    ctmp2->copy_block(0, offset + nact_ + j, n, 1, slice(offset +j*2+1, offset + j*2+2));
  }
  offset = (nclosed_+nact_)*2;
  // virtual (including positrons)
  for (int j=0; j!=nvirt_rel(); ++j) {
    ctmp2->copy_block(0, offset + j,            n, 1, slice(offset + j*2,   offset + j*2+1));
    ctmp2->copy_block(0, offset + nvirt_rel() + j,   n, 1, slice(offset + j*2+1, offset + j*2+2));
  }
  auto out = make_shared<RelCoeff_Block>(*ctmp2, nclosed_, nact_, nvirt_nr_, nneg_);
  return out;
}


std::shared_ptr<RelCoeff_Striped> RelCoeff_Block::striped_format() const {
  assert(nneg_ % 2 == 0);
  auto ctmp2 = clone();
  // Transforms a coefficient matrix from block format to striped format : assumes ordering is (c,a,v,positrons)
  // striped format
  int n = ndim();
  int offset = nclosed_;
  // closed
  for (int j=0; j!=nclosed_; ++j) {
    ctmp2->copy_block(0, j*2,   n, 1, slice(j, j+1));
    ctmp2->copy_block(0, j*2+1, n, 1, slice(offset + j, offset + j+1));
  }
  offset = nclosed_*2;
  // active
  for (int j=0; j!=nact_; ++j) {
    ctmp2->copy_block(0, offset + j*2,   n, 1, slice(offset + j,         offset + j+1));
    ctmp2->copy_block(0, offset + j*2+1, n, 1, slice(offset + nact_ + j, offset + nact_ + j+1));
  }
  offset = (nclosed_+nact_)*2;
  // vituals (including positrons)
  for (int j=0; j!=nvirt_rel(); ++j) {
    ctmp2->copy_block(0, offset + j*2,   n, 1, slice(offset + j,          offset + j+1));
    ctmp2->copy_block(0, offset + j*2+1, n, 1, slice(offset + nvirt_rel() + j, offset + nvirt_rel() + j+1));
  }
  auto out = make_shared<RelCoeff_Striped>(*ctmp2, nclosed_, nact_, nvirt_nr_, nneg_);
  return out;
}


std::shared_ptr<RelCoeff_Block> RelCoeff_Kramers::block_format() const {

  // TODO this could be optimized to avoid copying
  shared_ptr<ZMatrix> work = clone();
  array<shared_ptr<const ZMatrix>,2> tmp = {{ slice_copy(0, mdim()/2), slice_copy(mdim()/2, mdim()) }};

  int i = 0;
  work->copy_block(0, i, ndim(), nclosed_, tmp[0]->slice(0,nclosed_)); i += nclosed_;
  work->copy_block(0, i, ndim(), nclosed_, tmp[1]->slice(0,nclosed_)); i += nclosed_;
  work->copy_block(0, i, ndim(), nact_, tmp[0]->slice(nclosed_, nclosed_+nact_)); i += nact_;
  work->copy_block(0, i, ndim(), nact_, tmp[1]->slice(nclosed_, nclosed_+nact_)); i += nact_;
  work->copy_block(0, i, ndim(), nvirt_rel(), tmp[0]->slice(nclosed_+nact_, tmp[0]->mdim())); i += nvirt_rel();
  work->copy_block(0, i, ndim(), nvirt_rel(), tmp[1]->slice(nclosed_+nact_, tmp[1]->mdim()));

  auto out = make_shared<RelCoeff_Block>(*work, nclosed_, nact_, nvirt_nr_, nneg_);
  return out;
}

BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Striped)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Block)
