//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_coeff.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;


void ZCASSCF::init_kramers_coeff(shared_ptr<const ZMatrix> hcore, shared_ptr<const RelOverlap> overlap) {
  // Kramers-adapted coefficient via quaternion diagonalization
  const int nele = geom_->nele()-charge_;

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

  shared_ptr<ZMatrix> focktmp = make_shared<DFock>(geom_, hcore, coeff_->slice(0, nele), gaunt_, breit_, /*store_half*/false, /*robust*/false);
  quaternion(focktmp);

  shared_ptr<ZMatrix> s12 = overlap->tildex(1.0e-9);
  quaternion(s12);

  auto fock_tilde = make_shared<ZMatrix>(*s12 % (*focktmp) * *s12);


  // quaternion diagonalization
  {
    unique_ptr<double[]> eig(new double[fock_tilde->ndim()]);
    zquatev_(fock_tilde->ndim(), fock_tilde->data(), eig.get());
  }

  // re-order to kramers format and move negative energy states to virtual space
  auto ctmp = make_shared<ZMatrix>(*s12 * *fock_tilde);

  // rows: {L+, S+, L-, S-} -> {L+, L-, S+, S-}
  {
    assert(ctmp->ndim() % 4 == 0);
    const int n = ctmp->ndim()/4;
    const int m = ctmp->mdim();
    shared_ptr<ZMatrix> scratch = ctmp->get_submatrix(n, 0, n*2, m);
    ctmp->copy_block(n,   0, n, m, scratch->get_submatrix(n, 0, n, m));
    ctmp->copy_block(n*2, 0, n, m, scratch->get_submatrix(0, 0, n, m));
  }

  // move_positronic_orbitals;
  {
    auto move_one = [this, &ctmp](const int offset, const int block1, const int block2) {
      shared_ptr<ZMatrix> scratch = make_shared<ZMatrix>(ctmp->ndim(), block1+block2);
      scratch->copy_block(0,      0, ctmp->ndim(), block2, ctmp->slice(offset+block1, offset+block1+block2));
      scratch->copy_block(0, block2, ctmp->ndim(), block1, ctmp->slice(offset,        offset+block1));
      ctmp->copy_block(0, offset, ctmp->ndim(), block1+block2, scratch);
    };
    const int nneg2 = nneg_/2;
    move_one(           0, nneg2, nocc_+nvirt_-nneg2);
    move_one(nocc_+nvirt_, nneg2, nocc_+nvirt_-nneg2);
  }

  array<shared_ptr<const ZMatrix>,2> tmp = {{ ctmp->slice(0, ctmp->mdim()/2), ctmp->slice(ctmp->mdim()/2, ctmp->mdim()) }};

  shared_ptr<ZMatrix> ctmp2 = coeff_->clone();
  int i = 0;
  ctmp2->copy_block(0, i, coeff_->ndim(), nclosed_, tmp[0]->slice(0,nclosed_)); i += nclosed_;
  ctmp2->copy_block(0, i, coeff_->ndim(), nclosed_, tmp[1]->slice(0,nclosed_)); i += nclosed_;
  ctmp2->copy_block(0, i, coeff_->ndim(), nact_, tmp[0]->slice(nclosed_, nocc_)); i += nact_;
  ctmp2->copy_block(0, i, coeff_->ndim(), nact_, tmp[1]->slice(nclosed_, nocc_)); i += nact_;
  ctmp2->copy_block(0, i, coeff_->ndim(), nvirt_, tmp[0]->slice(nocc_, nocc_+nvirt_)); i += nvirt_;
  ctmp2->copy_block(0, i, coeff_->ndim(), nvirt_, tmp[1]->slice(nocc_, nocc_+nvirt_));
  coeff_ = ctmp2;
}


void ZCASSCF::kramers_adapt(shared_ptr<ZMatrix> o) const {
  // function to enforce time-reversal symmetry
  //    for a complex matrix o, that is SYMMETRIC under time reversal

  auto kramers_adapt_block = [this](shared_ptr<ZMatrix> o, unsigned int tfac, int nq1, int nq2, int joff, int ioff) {
    // function to enforce time-reversal symmetry for a given block of a complex matrix o.
    // tfac                 symmetry factor under time-reversal (symmetric : t=1, antisymm : t=-1)
    // nq1, nq2             nclosed_, nact_, nvirt_
    // ioff, joff           column and row offsets, respectively
    assert( o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt_)*2 == o->ndim() );
    assert( tfac == 1 || tfac == -1);
    const double t = tfac == 1 ? 1.0 : -1.0;
    for (int i = 0; i != nq2; ++i) {
      for (int j = 0; j != nq1; ++j) {
        // Diagonal contributions : "A" matrices in notation of T. Suae thesis
        o->element(joff+j, ioff+i) = ( o->element(joff+j, ioff+i) + conj(o->element(joff+j+nq1, ioff+i+nq2)) ) * 0.5;
        o->element(joff+j+nq1,ioff+i+nq2) = t * conj(o->element(joff+j,ioff+i));

        // Diagonal contributions : "B" matrices in notation of T. Suae thesis
        o->element(joff+nq1+j, ioff+i) = ( o->element(joff+nq1+j, ioff+i) - conj(o->element(joff+j, ioff+nq2+i)) ) * 0.5;
        o->element(joff+j, ioff+nq2+i) = -t * conj(o->element(joff+nq1+j,ioff+i));
      }
    }
  };

  assert(o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt_)*2 == o->ndim());
  const array<int,3> a0 {{nclosed_, nact_, nvirt_}};
  const array<int,3> a1 {{0, 2*nclosed_, 2*nocc_}};
  for (int ii = 0; ii !=3; ++ii) {
    for (int jj = 0; jj !=3; ++jj) {
      kramers_adapt_block(o,1,a0[jj],a0[ii],a1[jj],a1[ii]);
    }
  }
}

