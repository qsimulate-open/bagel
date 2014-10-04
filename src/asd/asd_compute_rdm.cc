//
// BAGEL - Parallel electron correlation program.
// Filename: asd_compute_rdm.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/asd/asd_base.h>
#include <src/asd/state_tensor.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;
using namespace btas;

//***************************************************************************************************************
void
ASD_base::compute_rdm () {
//***************************************************************************************************************
//const int norbA = dimer_->active_refs().first->nact();
//const int norbB = dimer_->active_refs().second->nact();

  // compute transformed gammas (J',J,zeta)
  StateTensor st(adiabats_, subspaces_base());
  st.print();

  // TODO parallelize
  // Loop over both tensors and mupltiply
  const int istate = 0;
  GammaTensor half;
  for (auto& i : *gammatensor_[0]) {
    for (auto& j : st) {
      // if the third index of the gamma tensor is identical to the first one of the state tensor we contract
      auto& ikey = i.first;
      auto& jkey = j.first;
      if (get<0>(jkey) == istate && get<2>(ikey) == get<1>(jkey)) {
        auto tag = make_tuple(get<0>(ikey), get<1>(ikey), get<2>(jkey));
        if (half.exist(tag)) {
          contract(1.0, *i.second, {0,1,2}, j.second, {1,3}, 0.0, *half.get_block(tag), {0,3,2});
        } else {
          auto data = make_shared<Tensor3<double>>(get<1>(ikey).nstates(), get<2>(jkey).nstates(), i.second->extent(2));
          contract(1.0, *i.second, {0,1,2}, j.second, {1,3}, 0.0, *data, {0,3,2});
          half.emplace(tag, data);
        }
      }
    }
  }

  worktensor_ = make_shared<GammaTensor>();
  for (auto& i : half) {
    for (auto& j : st) {
      auto& ikey = i.first;
      auto& jkey = j.first;
      if (get<0>(jkey) == istate && get<1>(ikey) == get<1>(jkey)) {
        auto tag = make_tuple(get<0>(ikey), get<2>(jkey), get<2>(ikey));

        // TODO check if this transformation is necessary
        if (worktensor_->exist(tag)) {
          contract(1.0, *i.second, {0,1,2}, j.second, {0,3}, 1.0, *worktensor_->get_block(tag), {3,1,2});
        } else {
          auto data = make_shared<Tensor3<double>>(get<2>(jkey).nstates(), get<2>(ikey).nstates(), i.second->extent(2));
          contract(1.0, *i.second, {0,1,2}, j.second, {0,3}, 0.0, *data, {3,1,2});
          worktensor_->emplace(tag, data);
        }
      }
    }
  }

/*
  //compute transformed gamma (I',I,nu)
  //first half
  GammaTensor half2;
  for (auto& i : *gammatensor_[1]) {
    for (auto& j : st) {
      auto& ikey = i.first;
      auto& jkey = j.first;
      if (get<0>(jkey) == istate &&  //ground state
          get<1>(ikey) == get<2>(jkey)) { //gamma (nu,[J'],J) == C(0,I',[J'])
        auto tag = make_tuple(get<0>(ikey), get<1>(jkey), get<2>(ikey)); // (nu,I',J)
        if (half2.exist(tag)) {
          contract(1.0, *i.second, {0,1,2}, j.second, {3,0}, 0.0, *half2.get_block(tag), {3,1,2});
        } else {
          auto data = make_shared<Tensor3<double>>(get<1>(jkey).nstates(), get<2>(ikey).nstates(), i.second->extent(2)); // (I',J,nu)
          contract(1.0, *i.second, {0,1,2}, j.second, {3,0}, 0.0, *data, {3,1,2});
          // Gamma(J',J,nu) * C(I',J') -> Gamma'(I',J,nu)
          //       0  1 2       3  0             3  1 2
          half2.emplace(tag, data);
        }
      }
    }
  }

  worktensor2_ = make_shared<GammaTensor>();
  for (auto& i : half2) {
    for (auto& j : st) {
      auto& ikey = i.first;
      auto& jkey = j.first;
      if (get<0>(jkey) == istate && //ground state
          get<2>(ikey) == get<2>(jkey)) { //gamma'(nu,I',[J]) == C(0,I,[J])
        auto tag = make_tuple(get<0>(ikey), get<1>(ikey), get<1>(jkey)); // (nu,I',I)

        // TODO check if this transformation is necessary
        if (worktensor2_->exist(tag)) {
          contract(1.0, *i.second, {0,1,2}, j.second, {3,1}, 1.0, *worktensor2_->get_block(tag), {0,3,2});
        } else {
          auto data = make_shared<Tensor3<double>>(get<1>(ikey).nstates(), get<1>(jkey).nstates(), i.second->extent(2)); // (I',I,nu)
          contract(1.0, *i.second, {0,1,2}, j.second, {3,1}, 0.0, *data, {0,3,2});
          // Gamma'(I',J,nu) * C(I,J) -> Gamma''(I',I,nu)
          //        0  1 2       3 1             0  3 2
          worktensor2_->emplace(tag, data);
        }
      }
    }
  }
*/

  const auto subspaces = subspaces_base();


  // diagonal subspaces
  for (auto& subspace : subspaces) {
    shared_ptr<RDM<1>> r1;
    shared_ptr<RDM<2>> r2;
    tie(r1,r2) = compute_diagonal_block_RDM(subspace);
    if (r1) assert(false); //*onerdm_ += *r1;
    if (r2) *twordm_ += *r2;
  }
  
  // off diagonal subspaces
  for (auto iAB = subspaces.begin(); iAB != subspaces.end(); ++iAB) {
    for (auto jAB = subspaces.begin(); jAB != iAB; ++jAB) {
      shared_ptr<RDM<1>> r1;
      shared_ptr<RDM<2>> r2;
      tie(r1,r2) = couple_blocks_RDM(*jAB, *iAB); //Lower-triangular (i<->j)
      if (r1) *onerdm_ += *r1;
      if (r2) *twordm_ += *r2;
    }
  }
  

  cout << "!@# Unsymmetrized 1RDM" << endl;
  onerdm_->print(1.0e-6);

  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();
  const int nactT = nactA + nactB;  

  //Symmetrize: D_AB (calculated) D_BA (uncalc.& symmetrized here)
  auto matBA = std::make_shared<Matrix>(nactB,nactA); //D_BA empty
  {
    auto low = {0, nactA};
    auto up  = {nactA, nactT};
    auto view = btas::make_view(onerdm_->range().slice(low,up), onerdm_->storage()); //D_AB sector of D (read ptr)
    auto matAB = std::make_shared<Matrix>(nactA,nactB); //D_AB empty
    std::copy(view.begin(), view.end(), matAB->begin()); //D_AB filled
    SMITH::sort_indices<1,0, 0,1, 1,1>(matAB->data(), matBA->data(), nactA, nactB); // transpose and fill D_BA
  }
  {
    auto low = {nactA, 0};
    auto up  = {nactT, nactA};
    auto outv = btas::make_rwview(onerdm_->range().slice(low,up), onerdm_->storage()); //D_BA sector of D (read & write ptr)
    std::copy(matBA->begin(), matBA->end(), outv.begin()); //copy D_BA -> D_BA sector of D
  }

  cout << "!@# Symmetrized 1RDM" << endl;
  onerdm_->print(1.0e-6);

  //Symmetrize: d(ABAA) note p18B
  {
    auto low = {0,nactA,0,0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactA); //empty d_ABAA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAA filled
    { //d(AAAB)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactA*nactB); //empty d_AAAB
      SMITH::sort_indices<0,3,2,1, 0,1, -1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AAAB
      auto low = {0,0,0,nactA};
      auto up  = {nactA,nactA,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_AAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AAAB into d_AAAB sector of d
    } 
    { //d(BAAA)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactA); //empty d_BAAA
      SMITH::sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_BAAA
      auto low = {nactA,0,0,0};
      auto up  = {nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAA into d_BAAA sector of d
    } 
    { //d(AABA)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactB*nactA); //empty d_AABA
      SMITH::sort_indices<3,0,1,2, 0,1, -1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AABA
      auto low = {0,0,nactA,0};
      auto up  = {nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_AABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AABA into d_AABA sector of d
    } 
  }
 
  //Symmetrize: d(ABBB) note p18B
  {
    auto low = {0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactB); //empty d_ABBB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBB filled
    { //d(BBAB)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactA*nactB); //empty d_BBAB
      SMITH::sort_indices<2,1,0,3, 0,1, -1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBAB
      auto low = {nactA,nactA,0,nactA};
      auto up  = {nactT,nactT,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAB into d_BBAB sector of d
    } 
    { //d(BABB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactB*nactB); //empty d_BABB
      SMITH::sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BABB
      auto low = {nactA,0,nactA,nactA};
      auto up  = {nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BABB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BABB into d_BABB sector of d
    } 
    { //d(BBBA)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactB*nactA); //empty d_BBBA
      SMITH::sort_indices<1,2,3,0, 0,1, -1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBBA
      auto low = {nactA,nactA,nactA,0};
      auto up  = {nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBBA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBBA into d_BBBA sector of d
    } 
  }


  //Symmetrize: d(ABBA) note p19
  {
    auto low = {0,nactA,nactA,0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactA); //empty d_ABBA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBA filled
    { //d(AABB)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactB*nactB); //empty d_AABB
      SMITH::sort_indices<0,3,2,1, 0,1, -1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactA); //reorder and fill d_AABB
      auto low = {0,0,nactA,nactA};
      auto up  = {nactA,nactA,nactT,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_AABB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AABB into d_AABB sector of d
    } 
    { //d(BBAA)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactA*nactA); //empty d_BBAA
      SMITH::sort_indices<2,1,0,3, 0,1, -1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactA); //reorder and fill d_BBAA
      auto low = {nactA,nactA,0,0};
      auto up  = {nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAA into d_BBAA sector of d
    } 
    { //d(BAAB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactB); //empty d_BAAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactA); //reorder and fill d_BAAB
      auto low = {nactA,0,0,nactA};
      auto up  = {nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAB into d_BAAB sector of d
    } 
  }

  //Symmetrize: d(ABAB) note p19
  {
    auto low = {0,nactA,0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactB); //empty d_ABAB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAB filled
    { //d(BABA)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactB*nactA); //empty d_BABA
      SMITH::sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactB); //reorder and fill d_BABA
      auto low = {nactA,0,nactA,0};
      auto up  = {nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBBA into d_BABA sector of d
    } 
  }


  //Energy calculation
  cout << "!@# Energy calculated from RDM:" << endl;
  const int nclosedA = dimer_->active_refs().first->nclosed();
  const int nclosedB = dimer_->active_refs().second->nclosed();
  cout << "Number of closed orbitals: A(" << nclosedA << "), B(" << nclosedB << ")" << endl;
  cout << "Number of active orbitals: A(" << nactA << "), B(" << nactB << ")" << endl;

  shared_ptr<const Matrix> ha = jop_->monomer_jop<0>()->mo1e()->matrix(); //h_AA
  shared_ptr<const Matrix> hb = jop_->monomer_jop<1>()->mo1e()->matrix(); //h_BB
  //                                                    CSymMatrix -> Matrix conversion
  shared_ptr<const Matrix> hc = jop_->cross_mo1e(); //h_AB

  auto int1 = make_shared<Matrix>(nactT,nactT);
  int1->zero();
  int1->copy_block(0,0,ha->ndim(),ha->mdim(),ha);
  int1->copy_block(nactA,nactA,hb->ndim(),hb->mdim(),hb);
  int1->copy_block(0,nactA,hc->ndim(),hc->mdim(),hc);
  int1->copy_block(nactA,0,hc->mdim(),hc->ndim(),hc->transpose());
  int1->print("1e integral",nactT);

  auto rdm1 = onerdm_->rdm1_mat(0);
  rdm1->print("1RDM",nactT);

}

//int n = 3;
//Matrix K(n,n);
//K.unit();
//btas::CRange<2> range(n*n, 1);
//const MatView gammaview(btas::make_view(range, K.storage()), /*localized*/true);
//
//auto MAT = make_shared<Matrix>( gammaview );
//MAT->print("a",n*n);
//cout << "(nxm)" << MAT->ndim() << MAT->mdim() << endl;

