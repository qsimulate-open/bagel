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

#include <src/multi/zcasscf/zcasscf.h>
#include <src/scf/dhf/dfock.h>
#include <src/util/math/quatmatrix.h>

using namespace std;
using namespace bagel;


// TODO This function uses QuatMatrix::diagonalize() - make sure the matrix is in the correct format
shared_ptr<RelCoeff_Striped> ZCASSCF::init_kramers_coeff_nonrel() {
  mute_stdcout();
  // Kramers-adapted coefficient via quaternion diagonalization
  const int nele = geom_->nele()-charge_;
  const int nocc = nclosed_ + nact_;
  const int nneg2 = coeff_->mdim()/4;
  const int nvirt_rel = coeff_->mdim()/2 - nocc;
  const int nvirt_nr = nvirt_rel - nneg2;

  shared_ptr<ZMatrix> coefftmp = nonrel_to_relcoeff(false);

  shared_ptr<ZMatrix> focktmp;
  shared_ptr<ZMatrix> ctmp;
  int norb = nele;//-2;
  ctmp = make_shared<ZMatrix>(coefftmp->ndim(), norb);
  ctmp->copy_block(0, 0, coefftmp->ndim(), norb/2, coefftmp->slice(0, norb/2));
  ctmp->copy_block(0, norb/2, coefftmp->ndim(), norb/2, coefftmp->slice(coefftmp->mdim()/2, coefftmp->mdim()/2+norb/2));
  focktmp = make_shared<DFock>(geom_, hcore_, ctmp, gaunt_, breit_, /*store_half*/false, /*robust*/false);

  shared_ptr<ZMatrix> fmo;
  if (tsymm_) {
    fmo = make_shared<QuatMatrix>(*coefftmp % (*focktmp) * *coefftmp);
#ifndef NDEBUG
    auto quatfmo = static_pointer_cast<const QuatMatrix>(fmo);
    const double tsymm_err = quatfmo->check_t_symmetry();
    if (tsymm_err > 1.0e-8)
      cout << "   ** Caution:  poor Kramers symmetry in fmo (ZCASSCF initialization) - error = " << scientific << setprecision(4) << tsymm_err << endl;
#endif
  } else {
    fmo = make_shared<ZMatrix>(*coefftmp % (*focktmp) * *coefftmp);
  }

  // quaternion diagonalization
  VectorB eig(fmo->ndim());
  fmo->diagonalize(eig);

  if (!tsymm_)
    RelCoeff::rearrange_eig(eig, fmo);

  auto ktmp = make_shared<RelCoeff_Kramers>(*fmo, nclosed_, nact_, nvirt_nr, nneg2*2);
  auto out = ktmp->block_format()->striped_format();
  resume_stdcout();
  return out;
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


void ZCASSCF::kramers_adapt(shared_ptr<ZRotFile> o, const int nclosed, const int nact, const int nvirt) {
  for (int i = 0; i != nclosed; ++i) {
    for (int j = 0; j != nvirt; ++j) {
      o->ele_vc(j, i) = (o->ele_vc(j, i) + conj(o->ele_vc(j+nvirt, i+nclosed))) * 0.5;
      o->ele_vc(j+nvirt, i+nclosed) = conj(o->ele_vc(j, i));

      o->ele_vc(j+nvirt, i) = (o->ele_vc(j+nvirt, i) - conj(o->ele_vc(j, i+nclosed))) * 0.5;
      o->ele_vc(j, i+nclosed) = - conj(o->ele_vc(j+nvirt, i));
    }
  }
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nvirt; ++j) {
      o->ele_va(j, i) = (o->ele_va(j, i) + conj(o->ele_va(j+nvirt, i+nact))) * 0.5;
      o->ele_va(j+nvirt, i+nact) = conj(o->ele_va(j, i));

      o->ele_va(j+nvirt, i) = (o->ele_va(j+nvirt, i) - conj(o->ele_va(j, i+nact))) * 0.5;
      o->ele_va(j, i+nact) = - conj(o->ele_va(j+nvirt, i));
    }
  }
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nclosed; ++j) {
      o->ele_ca(j, i) = (o->ele_ca(j, i) + conj(o->ele_ca(j+nclosed, i+nact))) * 0.5;
      o->ele_ca(j+nclosed, i+nact) = conj(o->ele_ca(j, i));

      o->ele_ca(j+nclosed, i) = (o->ele_ca(j+nclosed, i) - conj(o->ele_ca(j, i+nact))) * 0.5;
      o->ele_ca(j, i+nact) = - conj(o->ele_ca(j+nclosed, i));
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
  VectorB eig2(shalf->mdim());
  shalf->diagonalize(eig2);
  for (int k = 0; k != shalf->mdim(); ++k) {
    if (real(eig2[k]) >= 0.0)
      blas::scale_n(sqrt(sqrt(eig2(k))), shalf->element_ptr(0, k), shalf->ndim());
  }
  *shalf = *shalf ^ *shalf;
  auto tcoeff = make_shared<ZMatrix>(nr_coeff_->ndim(), nr_coeff_->mdim());
  tcoeff->add_real_block(1.0, 0, 0, n, n, *nr_coeff_);
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
