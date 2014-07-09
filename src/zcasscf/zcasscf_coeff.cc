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


void ZCASSCF::init_kramers_coeff() {
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

  shared_ptr<ZMatrix> coefftmp;
  if (nr_coeff_ != nullptr)
    coefftmp = nonrel_to_relcoeff(false);

  shared_ptr<ZMatrix> focktmp;
  shared_ptr<ZMatrix> ctmp;
  if (nr_coeff_ == nullptr) {
    int norb = nele;// - 2 > 0 ? nele-2 : nele;
    focktmp = make_shared<DFock>(geom_, hcore_, coeff_->slice_copy(0, norb), gaunt_, breit_, /*store_half*/false, /*robust*/false);
    quaternion(focktmp);

    shared_ptr<ZMatrix> s12 = overlap_->tildex(1.0e-9);
    quaternion(s12);

    auto fock_tilde = make_shared<ZMatrix>(*s12 % (*focktmp) * *s12);

    // quaternion diagonalization
    {
      unique_ptr<double[]> eig(new double[fock_tilde->ndim()]);
      zquatev_(fock_tilde->ndim(), fock_tilde->data(), eig.get());
    }

    // re-order to kramers format and move negative energy states to virtual space
    ctmp = make_shared<ZMatrix>(*s12 * *fock_tilde);

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
  } else {//if (nele%2 == 0 && nele - 2 > 0) {
    int norb = nele;//-2;
    auto ctmp = make_shared<ZMatrix>(coefftmp->ndim(), norb);
    ctmp->copy_block(0, 0, coefftmp->ndim(), norb/2, coefftmp->slice(0, norb/2));
    ctmp->copy_block(0, norb/2, coefftmp->ndim(), norb/2, coefftmp->slice(coefftmp->mdim()/2, coefftmp->mdim()/2+norb/2));
    focktmp = make_shared<DFock>(geom_, hcore_, ctmp, gaunt_, breit_, /*store_half*/false, /*robust*/false);
    auto fmo = make_shared<ZMatrix>(*coefftmp % *focktmp * *coefftmp);
    // quaternion diagonalization
    {
      unique_ptr<double[]> eig(new double[fmo->ndim()]);
      zquatev_(fmo->ndim(), fmo->data(), eig.get());
      // move_positronic_orbitals;
      {
        auto move_one = [this, &fmo](const int offset, const int block1, const int block2) {
          shared_ptr<ZMatrix> scratch = make_shared<ZMatrix>(fmo->ndim(), block1+block2);
          scratch->copy_block(0,      0, fmo->ndim(), block2, fmo->slice(offset+block1, offset+block1+block2));
          scratch->copy_block(0, block2, fmo->ndim(), block1, fmo->slice(offset,        offset+block1));
          fmo->copy_block(0, offset, fmo->ndim(), block1+block2, scratch);
        };
        const int nneg2 = nneg_/2;
        move_one(           0, nneg2, nocc_+nvirt_-nneg2);
        move_one(nocc_+nvirt_, nneg2, nocc_+nvirt_-nneg2);
      }
      coefftmp = make_shared<ZMatrix>(*coefftmp * *fmo);
    }
  }

  array<shared_ptr<const ZMatrix>,2> tmp;
  if (nr_coeff_ == nullptr) {
    tmp = {{ ctmp->slice_copy(0, ctmp->mdim()/2), ctmp->slice_copy(ctmp->mdim()/2, ctmp->mdim()) }};
  } else{
    tmp = {{ coefftmp->slice_copy(0, coefftmp->mdim()/2), coefftmp->slice_copy(coefftmp->mdim()/2, coefftmp->mdim()) }};
  }
  shared_ptr<ZMatrix> ctmp2 = coeff_->clone();

  { // striped format
    int n = coeff_->ndim();
    // closed
    for (int j=0; j!=nclosed_; ++j) {
      ctmp2->copy_block(0, j*2,   n, 1, tmp[0]->slice(j, j+1));
      ctmp2->copy_block(0, j*2+1, n, 1, tmp[1]->slice(j, j+1));
    }
    int offset = nclosed_*2;
    // active
    for (int j=0; j!=nact_; ++j) {
      ctmp2->copy_block(0, offset + j*2,   n, 1, tmp[0]->slice(offset/2 + j, offset/2 + j+1));
      ctmp2->copy_block(0, offset + j*2+1, n, 1, tmp[1]->slice(offset/2 + j, offset/2 + j+1));
    }
    offset = nocc_*2;
    // active
    for (int j=0; j!=nvirtnr_; ++j) {
      ctmp2->copy_block(0, offset + j*2,   n, 1, tmp[0]->slice(offset/2 + j, offset/2 + j+1));
      ctmp2->copy_block(0, offset + j*2+1, n, 1, tmp[1]->slice(offset/2 + j, offset/2 + j+1));
    }
    offset = ctmp2->mdim()/2;
    // positrons
    for (int j=0; j!=nneg_/2; ++j) {
      ctmp2->copy_block(0, offset + j*2,   n, 1, tmp[0]->slice(nocc_+nvirtnr_ + j, nocc_+nvirtnr_ + j+1));
      ctmp2->copy_block(0, offset + j*2+1, n, 1, tmp[1]->slice(nocc_+nvirtnr_ + j, nocc_+nvirtnr_ + j+1));
    }
  }
  coeff_ = ctmp2;
}


void ZCASSCF::kramers_adapt(shared_ptr<ZMatrix> o, const int nvirt) const {
  // function to enforce time-reversal symmetry
  //    for a complex matrix o, that is SYMMETRIC under time reversal

  auto kramers_adapt_block = [this](shared_ptr<ZMatrix> o, unsigned int tfac, int nq1, int nq2, int joff, int ioff, int nvirt) {
    // function to enforce time-reversal symmetry for a given block of a complex matrix o.
    // tfac                 symmetry factor under time-reversal (symmetric : t=1, antisymm : t=-1)
    // nq1, nq2             nclosed_, nact_, nvirt
    // ioff, joff           column and row offsets, respectively
    assert( o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt)*2 == o->ndim() );
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

  assert(o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt)*2 == o->ndim());
  const array<int,3> a0 {{nclosed_, nact_, nvirt}};
  const array<int,3> a1 {{0, 2*nclosed_, 2*nocc_}};
  for (int ii = 0; ii !=3; ++ii) {
    for (int jj = 0; jj !=3; ++jj) {
      kramers_adapt_block(o,1,a0[jj],a0[ii],a1[jj],a1[ii],nvirt);
    }
  }
}


void ZCASSCF::kramers_adapt(shared_ptr<ZRotFile> o, const int nvirt) const {
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt; ++j) {
      o->ele_vc(j, i) = (o->ele_vc(j, i) + conj(o->ele_vc(j+nvirt, i+nclosed_))) * 0.5;
      o->ele_vc(j+nvirt, i+nclosed_) = conj(o->ele_vc(j, i));

      o->ele_vc(j+nvirt, i) = (o->ele_vc(j+nvirt, i) - conj(o->ele_vc(j, i+nclosed_))) * 0.5;
      o->ele_vc(j, i+nclosed_) = - conj(o->ele_vc(j+nvirt, i));
    }
  }
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt; ++j) {
      o->ele_va(j, i) = (o->ele_va(j, i) + conj(o->ele_va(j+nvirt, i+nact_))) * 0.5;
      o->ele_va(j+nvirt, i+nact_) = conj(o->ele_va(j, i));

      o->ele_va(j+nvirt, i) = (o->ele_va(j+nvirt, i) - conj(o->ele_va(j, i+nact_))) * 0.5;
      o->ele_va(j, i+nact_) = - conj(o->ele_va(j+nvirt, i));
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


shared_ptr<ZMatrix> ZCASSCF::nonrel_to_relcoeff(const bool stripes) const {
  // constructs a relativistic coefficient for electronic componenets from a non-rel coefficient
  int n = nr_coeff_->ndim();
  int nvirtnr = nvirt_ - nneg_/2;
  auto ctmp = make_shared<ZMatrix>(n*4, n*4);

  // compute T^(-1/2) and S^(1/2)
  auto t12 = overlap_->get_submatrix(n*2, n*2, n, n)->copy();
  t12->inverse_half(1.0e-10);
  auto shalf = overlap_->get_submatrix(0, 0, n, n)->copy();
  unique_ptr<double[]> eig2(new double[shalf->mdim()]);
  shalf->diagonalize(eig2.get());
  for (int k = 0; k != shalf->mdim(); ++k) {
    if (real(eig2[k]) >= 0.0)
      blas::scale_n(sqrt(sqrt(eig2[k])), shalf->element_ptr(0, k), shalf->ndim());
  }
  *shalf = *shalf ^ *shalf;
  auto tcoeff = make_shared<ZMatrix>(nr_coeff_->ndim(), nr_coeff_->mdim());
  tcoeff->add_real_block(1.0, 0, 0, n, n, nr_coeff_->data());
  *tcoeff = *t12 * *shalf * *tcoeff;
  // copy non rel values into "+/-" stripes
  if (stripes) {
    // closed
    for (int j=0; j!=nclosed_; ++j) {
      ctmp->add_real_block(1.0, 0, j*2,   n, 1, nr_coeff_->slice(j, j+1));
      ctmp->add_real_block(1.0, n, j*2+1, n, 1, nr_coeff_->slice(j, j+1));
    }
    // active
    for (int j=0; j!=nact_; ++j) {
      ctmp->add_real_block(1.0, 0, nclosed_*2 + j*2,   n, 1, nr_coeff_->slice(nclosed_ + j, nclosed_ + j+1));
      ctmp->add_real_block(1.0, n, nclosed_*2 + j*2+1, n, 1, nr_coeff_->slice(nclosed_ + j, nclosed_ + j+1));
    }
    // virtuals
    for (int j=0; j!=nvirtnr; ++j) {
      ctmp->add_real_block(1.0, 0, nocc_*2 + j*2,   n, 1, nr_coeff_->slice(nocc_ + j, nocc_ + j+1));
      ctmp->add_real_block(1.0, n, nocc_*2 + j*2+1, n, 1, nr_coeff_->slice(nocc_ + j, nocc_ + j+1));
    }
    // positrons
    for (int j=0; j!=nneg_/2; ++j) {
      ctmp->copy_block(n*2, ctmp->mdim() - j*2-2,   n, 1, tcoeff->slice(j, j+1));
      ctmp->copy_block(n*3, ctmp->mdim() - j*2-1, n, 1, tcoeff->slice(j, j+1));
    }
  } else { // block structure
    int m2 = ctmp->mdim()/2;
    // closed
    for (int j=0; j!=nclosed_; ++j) {
      ctmp->add_real_block(1.0, 0, j,   n, 1, nr_coeff_->slice(j, j+1));
      ctmp->add_real_block(1.0, n, m2 + j, n, 1, nr_coeff_->slice(j, j+1));
    }
    // active
    for (int j=0; j!=nact_; ++j) {
      ctmp->add_real_block(1.0, 0, nclosed_ + j, n, 1, nr_coeff_->slice(nclosed_ + j, nclosed_ + j+1));
      ctmp->add_real_block(1.0, n, m2 + nclosed_ + j, n, 1, nr_coeff_->slice(nclosed_ + j, nclosed_ + j+1));
    }
    // virtuals
    for (int j=0; j!=nvirtnr; ++j) {
      ctmp->add_real_block(1.0, 0, nocc_ + j,   n, 1, nr_coeff_->slice(nocc_ + j, nocc_ + j+1));
      ctmp->add_real_block(1.0, n, m2 + nocc_ + j, n, 1, nr_coeff_->slice(nocc_ + j, nocc_ + j+1));
    }
    // positrons
    for (int j=0; j!=nneg_/2; ++j) {
      ctmp->copy_block(n*2, m2 - j-1,   n, 1, tcoeff->slice(j, j+1));
      ctmp->copy_block(n*3, ctmp->mdim() - j-1, n, 1, tcoeff->slice(j, j+1));
    }
  }
  return ctmp;
}


// Transforms a coefficient matrix from striped format to block format : assumes ordering is (c,a,v,positrons)
shared_ptr<ZMatrix> ZCASSCF::coeff_stripe_to_block(shared_ptr<const ZMatrix> coeff) const {
  assert(coeff->ndim() == coeff->mdim());
  auto ctmp2 = coeff->clone();
  { // block format
    int n = coeff->ndim();
    // closed
    for (int j=0; j!=nclosed_; ++j) {
      ctmp2->copy_block(0,            j, n, 1, coeff->slice(j*2  , j*2+1));
      ctmp2->copy_block(0, nclosed_ + j, n, 1, coeff->slice(j*2+1, j*2+2));
    }
    int offset = nclosed_*2;
    // active
    for (int j=0; j!=nact_; ++j) {
      ctmp2->copy_block(0, offset + j,         n, 1, coeff->slice(offset +j*2,   offset + j*2+1));
      ctmp2->copy_block(0, offset + nact_ + j, n, 1, coeff->slice(offset +j*2+1, offset + j*2+2));
    }
    offset = nocc_*2;
    // virtual (including positrons)
    for (int j=0; j!=nvirt_; ++j) {
      ctmp2->copy_block(0, offset + j,            n, 1, coeff->slice(offset + j*2,   offset + j*2+1));
      ctmp2->copy_block(0, offset + nvirt_ + j,   n, 1, coeff->slice(offset + j*2+1, offset + j*2+2));
    }
  }
   return ctmp2;
}


// Eliminates the positronic entries for the given rot file
void ZCASSCF::zero_positronic_elements(shared_ptr<ZRotFile> rot) {
  int nr_nvirt = nvirt_ - nneg_/2;
  for (int i = 0; i != nclosed_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_vc(j + nr_nvirt, i) =  complex<double> (0.0,0.0);
      rot->ele_vc(j + nr_nvirt + nvirt_, i) =  complex<double> (0.0,0.0);
    }
  }
  for (int i = 0; i != nact_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_va(j + nr_nvirt, i) =  complex<double> (0.0,0.0);
      rot->ele_va(j + nr_nvirt + nvirt_, i) =  complex<double> (0.0,0.0);
    }
  }
}
