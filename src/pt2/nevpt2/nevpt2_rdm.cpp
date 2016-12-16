//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: nevpt2_rdm.cc
// Copyright (C) 2014 Toru Shiozaki
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

#ifdef NEVPT2IMPL

template<>
void NEVPT2<double>::compute_rdm() {
  // rdm 1
  {
    auto tmp = ref_->rdm1(istate_)->rdm1_mat(/*nclosed_*/0, false);
    tmp->localize();
    rdm1_ = tmp;
  }
  // rdm 2
  {
    auto tmp = make_shared<MatType>(nact_*nact_, nact_*nact_, true);
    shared_ptr<const RDM<2>> r2 = ref_->rdm2(istate_);
    sort_indices<0,2,1,3,0,1,1,1>(r2->data(), tmp->data(), nact_, nact_, nact_, nact_);
    rdm2_ = tmp;
  }
  // rdm 3 and 4
  {
    auto tmp3 = make_shared<MatType>(nact_*nact_*nact_, nact_*nact_*nact_, true);
    auto tmp4 = make_shared<MatType>(nact_*nact_*nact_*nact_, nact_*nact_*nact_*nact_, true);
    shared_ptr<const RDM<3>> r3;
    shared_ptr<const RDM<4>> r4;
    tie(r3, r4) = ref_->rdm34(istate_, istate_);
    sort_indices<0,2,4,  1,3,5,  0,1,1,1>(r3->data(), tmp3->data(), nact_, nact_, nact_, nact_, nact_, nact_);
    sort_indices<0,2,4,6,1,3,5,7,0,1,1,1>(r4->data(), tmp4->data(), nact_, nact_, nact_, nact_, nact_, nact_, nact_, nact_);
    rdm3_ = tmp3;
    rdm4_ = tmp4;
  }
}

template<>
void NEVPT2<complex<double>>::compute_rdm() {
  auto ref = dynamic_pointer_cast<const RelReference>(ref_);
  // rdm 1
  {
    auto tmp = make_shared<MatType>(nact_, nact_, true);
    auto rdm1k = ref->rdm1(istate_, istate_);
    auto r1 = expand_kramers(rdm1k, nact_/2);
    copy_n(r1->data(), nact_*nact_, tmp->data());
    rdm1_ = tmp;
  }
  // rdm 2
  {
    auto tmp = make_shared<MatType>(nact_*nact_, nact_*nact_, true);
    auto rdm2k = ref->rdm2(istate_, istate_);
    auto r2 = expand_kramers(rdm2k, nact_/2);
    sort_indices<0,2,1,3,0,1,1,1>(r2->data(), tmp->data(), nact_, nact_, nact_, nact_);
    rdm2_ = tmp;
  }
  // rdm 3
  {
    auto tmp3 = make_shared<MatType>(nact_*nact_*nact_, nact_*nact_*nact_, true);
    auto rdm3k = ref->rdm3(istate_, istate_);
    auto r3 = expand_kramers(rdm3k, nact_/2);
    sort_indices<0,2,4,1,3,5,0,1,1,1>(r3->data(), tmp3->data(), nact_, nact_, nact_, nact_, nact_, nact_);
    rdm3_ = tmp3;
  }
  // TODO rdm 4 is too large to do this way - implement direct computation in ARDM3 later
  // rdm 4
  {
    auto tmp4 = make_shared<MatType>(nact_*nact_*nact_*nact_, nact_*nact_*nact_*nact_, true);
    auto rdm4k = ref->rdm4(istate_, istate_);
    auto r4 = expand_kramers(rdm4k, nact_/2);
    sort_indices<0,2,4,6,1,3,5,7,0,1,1,1>(r4->data(), tmp4->data(), nact_, nact_, nact_, nact_, nact_, nact_, nact_, nact_);
    rdm4_ = tmp4;
  }
}


template<typename DataType>
void NEVPT2<DataType>::compute_asrdm() {
  assert(rdm1_ && rdm2_ && rdm3_ && rdm4_);
  auto id2 = [this](                          const int k, const int l) { return         (        (k+nact_*l)); };
  auto id3 = [this](             const int j, const int k, const int l) { return         (j+nact_*(k+nact_*l)); };
  auto id4 = [this](const int i, const int j, const int k, const int l) { return i+nact_*(j+nact_*(k+nact_*l)); };

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;

  // amat = <a+ a b+ b>, <a+ a b+ b c+ c>, and <a+ a b+ b c+ c d+ d>
  // also srdm2 = <0|a+p bp cq d+q|0>
  shared_ptr<MatType> ardm2 = rdm2_->clone();
  shared_ptr<MatType> srdm2 = rdm2_->clone();
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j) {
      for (int k = 0; k != nact_; ++k)
        blas::ax_plus_y_n(1.0, rdm2_->element_ptr(nact_*j,k+nact_*i), nact_, ardm2->element_ptr(nact_*k,j+nact_*i));
      blas::ax_plus_y_n(1.0, rdm1_->element_ptr(0,i), nact_, ardm2->element_ptr(nact_*j,j+nact_*i));
      blas::ax_plus_y_n(fac2, rdm1_->element_ptr(0,j), nact_, srdm2->element_ptr(nact_*j,i+nact_*i));
    }
  sort_indices<0,2,1,1,1,-1,1>(ardm2->data(), srdm2->data(), nact_*nact_, nact_, nact_);
  shared_ptr<MatType> ardm3 = rdm3_->clone();
  shared_ptr<MatType> srdm3 = rdm3_->clone(); // <a+ a b b+ c+ c>
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k) {
        for (int l = 0; l != nact_; ++l) {
          for (int m = 0; m != nact_; ++m)
            blas::ax_plus_y_n(1.0,  rdm3_->element_ptr(id3(0,l,j),id3(m,k,i)), nact_, ardm3->element_ptr(id3(0,m,l),id3(k,j,i)));
          blas::ax_plus_y_n(1.0, ardm2->element_ptr(nact_*k,j+nact_*i), nact_, ardm3->element_ptr(id3(0,l,l),id3(k,j,i)));
          blas::ax_plus_y_n(1.0, rdm2_->element_ptr(nact_*k,l+nact_*i), nact_, ardm3->element_ptr(id3(0,l,k),id3(j,j,i)));
          blas::ax_plus_y_n(1.0, rdm2_->element_ptr(nact_*k,i+nact_*j), nact_, ardm3->element_ptr(id3(0,l,k),id3(j,l,i)));
        }
        blas::ax_plus_y_n(fac2, ardm2->element_ptr(0, id2(j,i)), nact_*nact_, srdm3->element_ptr(id3(0,0,k),id3(k,j,i)));
      }
  sort_indices<0,2,1,3,1,1,-1,1>(ardm3->data(), srdm3->data(), nact_*nact_, nact_, nact_, nact_*nact_);
  shared_ptr<MatType> ardm4 = rdm4_->clone();
  for (int h = 0; h != nact_; ++h)
    for (int g = 0; g != nact_; ++g)
      for (int f = 0; f != nact_; ++f)
        for (int e = 0; e != nact_; ++e)
          for (int d = 0; d != nact_; ++d) {
            blas::ax_plus_y_n(-1.0, ardm2->element_ptr(id2(0,f),id2(g,h)), nact_, ardm4->element_ptr(id4(0,d,d,e),id4(e,f,g,h)));
            blas::ax_plus_y_n(-1.0, ardm2->element_ptr(id2(0,d),id2(g,h)), nact_, ardm4->element_ptr(id4(0,e,f,d),id4(e,f,g,h)));
            for (int c = 0; c != nact_; ++c) {
              blas::ax_plus_y_n(1.0, ardm3->element_ptr(id3(0,d,e),id3(f,g,h)), nact_, ardm4->element_ptr(id4(0,c,c,d),id4(e,f,g,h)));
              blas::ax_plus_y_n(1.0, ardm3->element_ptr(id3(0,c,d),id3(f,g,h)), nact_, ardm4->element_ptr(id4(0,c,d,e),id4(e,f,g,h)));
              blas::ax_plus_y_n(1.0, ardm3->element_ptr(id3(0,f,c),id3(d,g,h)), nact_, ardm4->element_ptr(id4(0,e,c,d),id4(e,f,g,h)));
              blas::ax_plus_y_n(1.0, rdm3_->element_ptr(id3(0,c,e),id3(f,d,h)), nact_, ardm4->element_ptr(id4(0,f,c,d),id4(e,g,g,h)));
              blas::ax_plus_y_n(1.0, rdm3_->element_ptr(id3(0,c,e),id3(d,h,f)), nact_, ardm4->element_ptr(id4(0,d,c,g),id4(e,f,g,h)));
              blas::ax_plus_y_n(1.0, rdm3_->element_ptr(id3(0,c,e),id3(h,d,f)), nact_, ardm4->element_ptr(id4(0,g,c,d),id4(e,f,g,h)));
              for (int b = 0; b != nact_; ++b)
                blas::ax_plus_y_n(1.0, rdm4_->element_ptr(id4(0,c,e,g),id4(b,d,f,h)), nact_, ardm4->element_ptr(id4(0,b,c,d),id4(e,f,g,h)));
            }
          }
  ardm2_ = ardm2;
  ardm3_ = ardm3;
  ardm4_ = ardm4;
  srdm2_ = srdm2;
  srdm3_ = srdm3;
}


template<typename DataType>
void NEVPT2<DataType>::compute_hrdm() {
  assert(rdm1_ && rdm2_ && rdm3_ && srdm2_);

  auto id3 = [this](const int j, const int k, const int l) { return j+nact_*(k+nact_*l); };

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;

  shared_ptr<MatType> unit = rdm1_->clone(); unit->unit();
  shared_ptr<MatType> hrdm2 = rdm2_->copy();
  shared_ptr<const MatType> hrdm1 = make_shared<MatType>(*unit*fac2 - *rdm1_);

  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nact_; ++j) {
      blas::ax_plus_y_n(fac2, hrdm1->element_ptr(0,i), nact_, hrdm2->element_ptr(nact_*j, i+nact_*j));
      blas::ax_plus_y_n(1.0, rdm1_->element_ptr(0,i), nact_, hrdm2->element_ptr(nact_*j, j+nact_*i));
      for (int k = 0; k != nact_; ++k) {
        hrdm2->element(k+nact_*j, i+nact_*k) -= hrdm1->element(j,i);
        hrdm2->element(k+nact_*j, k+nact_*i) -= fac2 * rdm1_->element(j,i);
      }
    }
  }
  auto hrdm3 = rdm3_->clone();
  auto hrdm3tmp = rdm3_->clone();
  auto srdm2trans = srdm2_->transpose_conjg();
  auto rdm2reo = rdm2_->clone();
  sort_indices<1,0,2,0,1,1,1>(rdm2_->data(), rdm2reo->data(), nact_, nact_, nact_*nact_);
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k) {
        for (int m = 0; m != nact_; ++m) {
          blas::ax_plus_y_n(-1.0, hrdm2->element_ptr(nact_*k, i+nact_*j), nact_, hrdm3tmp->element_ptr(id3(0,m,k),id3(j,i,m)));
          blas::ax_plus_y_n(1.0, srdm2trans->element_ptr(nact_*j, i+nact_*k), nact_, hrdm3tmp->element_ptr(id3(0,m,k),id3(j,m,i)));
          blas::ax_plus_y_n(fac2, rdm2reo->element_ptr(nact_*k, i+nact_*j), nact_, hrdm3tmp->element_ptr(id3(0,m,k),id3(m,j,i)));
        }
      }
  sort_indices<1,0,2,0,1,1,1>(hrdm3tmp->data(), hrdm3->data(), nact_, nact_, nact_*nact_*nact_*nact_);

  hrdm3->ax_plus_y(-1.0, rdm3_);
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k) {
        blas::ax_plus_y_n(fac2, hrdm2->element_ptr(0, j+nact_*i), nact_*nact_, hrdm3->element_ptr(id3(0,0,k),id3(j,i,k)));
        blas::ax_plus_y_n(-1.0, rdm2_->element_ptr(0, i+nact_*j), nact_*nact_, hrdm3->element_ptr(id3(0,0,k),id3(k,j,i)));
        for (int l = 0; l != nact_; ++l) {
          blas::ax_plus_y_n(-1.0, hrdm2->element_ptr(nact_*k,j+nact_*i), nact_, hrdm3->element_ptr(id3(0,l,k),id3(j,i,l)));
          blas::ax_plus_y_n(1.0, srdm2trans->element_ptr(nact_*j, i+nact_*k), nact_, hrdm3->element_ptr(id3(0,k,l),id3(j,l,i)));
          blas::ax_plus_y_n(-fac2, srdm2trans->element_ptr(nact_*j, i+nact_*k), nact_, hrdm3->element_ptr(id3(0,l,k),id3(j,l,i)));
          blas::ax_plus_y_n(-1.0, rdm2reo->element_ptr(nact_*k, i+nact_*j), nact_, hrdm3->element_ptr(id3(0,l,k),id3(l,j,i)));
        }
      }
  hrdm1_ = hrdm1;
  assert(hrdm1_->is_hermitian());
  hrdm2_ = hrdm2;
  assert(hrdm2_->is_hermitian());
  hrdm3_ = hrdm3;
  assert(hrdm3_->is_hermitian());
}

#endif
