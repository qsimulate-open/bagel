//
// BAGEL - Parallel electron correlation program.
// Filename: asd_rdm3.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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
#include <src/smith/prim_op.h>
#include <cmath>

using namespace std;
using namespace bagel;

void ASD_base::initialize_3RDM() {
  cout << "Initialize 4RDM" << endl;
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();

  //E_ai,bj,ck,dl = sum c'b'a'ijk
  //#of B indices = 0
  rdm3_.emplace(string("monomerA"), make_shared<Matrix>(pow(nactA,6), 1)); //monomer A
  //# = 1
  rdm3_.emplace(string("k"), make_shared<Matrix>(pow(nactA,5), nactB)); // k
  //# = 2
  rdm3_.emplace(string("jk"), make_shared<Matrix>(pow(nactA,4), pow(nactB,2))); // jk
  rdm3_.emplace(string("ai"), make_shared<Matrix>(pow(nactA,4), pow(nactB,2))); // ai
  rdm3_.emplace(string("aj"), make_shared<Matrix>(pow(nactA,4), pow(nactB,2))); // aj
  //# = 3
  rdm3_.emplace(string("ijk"), make_shared<Matrix>(pow(nactA,3), pow(nactB,3))); // ijk
  rdm3_.emplace(string("aij"), make_shared<Matrix>(pow(nactA,3), pow(nactB,3))); // aij
  rdm3_.emplace(string("ajk"), make_shared<Matrix>(pow(nactA,3), pow(nactB,3))); // ajk
  //# = 4 
  rdm3_.emplace(string("aijk"), make_shared<Matrix>(pow(nactA,2), pow(nactB,4))); // aijk
  rdm3_.emplace(string("baij"), make_shared<Matrix>(pow(nactA,2), pow(nactB,4))); // baij
  rdm3_.emplace(string("bajk"), make_shared<Matrix>(pow(nactA,2), pow(nactB,4))); // bajk
  //# = 5
  rdm3_.emplace(string("baijk"), make_shared<Matrix>(nactA, pow(nactB,5))); //baijk
  //# = 6
  rdm3_.emplace(string("monomerB"), make_shared<Matrix>(pow(nactA,6), 1)); //monomer B

  cout << "# of nonredundnat Dimer 3RDM = " << rdm3_.size() << endl;

}
#if 0
double ASD_base::element_3RDM(const int a, const int i, const int b, const int j,const int c, const int k) {
  const int nactA = dimer_->embedded_refs().first->nact();
  vector<int>  left{c,b,a};
  vector<int>  right{k,j,i};
  vector<bool> bar_left{false,false,false};
  vector<bool> bar_right{false,false,false};
  int nB_left{0};
  //left
  for (int i=0; i<3; i++) {
    int temp = index - NactA;
    if(temp >= 0) { //B index
      left.at(i) = temp;
      bar_left.at(i) = true;
      ++nB_left;
    }
  }
  //right
  int nB_right{0};
  for (int i=0; i<3; i++) {
    int temp = index - NactA;
    if(temp >= 0) { //B index
      right.at(i) = temp;
      bar_right.at(i) = temp;
      ++nB_right;
    }
  }
  //swap
  if(nB_left > nB_right) {
    swap(left,right);
    swap(nB_left, nB_right);
    swap(bar_left, bar_right);
  } 
  //move
  for (int i=1; i<3; i++) { //skip the first element
    if(!bar_right.at(i)) continue; //A index 
    //move down B index
    assert(bar_right.at(i-1) == true); //error if i-1 = true
    for (int j=i; j=0; j--) {
    }
  }

  return ...
}
#endif

tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_diagonal_block_RDM34(const DimerSubspace_base& subspace) const {
  array<MonomerKey,4> keys {{ subspace.monomerkey<0>(), subspace.monomerkey<1>(), subspace.monomerkey<0>(), subspace.monomerkey<1>() }};
  auto out = compute_diag_RDM34(keys, /*subspace diagonal*/true);
  return out;
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>>
ASD_base::couple_blocks_RDM34(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {
//***************************************************************************************************************
  cout << "couple_block_RDM34 enetered.." << endl;
  Coupling term_type = coupling_type_RDM34(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);
  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }
  
  tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> out;
  std::array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = make_tuple(nullptr,nullptr); break;
    case Coupling::diagonal :
      out = compute_diag_RDM34(keys, /*subspace diagonal*/false); break;
    case Coupling::aET :
      out = compute_aET_RDM34(keys); break;
    case Coupling::bET :
      out = compute_bET_RDM34(keys); break;
    case Coupling::abFlip :
      out = compute_abFlip_RDM34(keys); break;
    case Coupling::abET :
      out = compute_abET_RDM34(keys); break;
    case Coupling::aaET :
      out = compute_aaET_RDM34(keys); break;
    case Coupling::bbET :
      out = compute_bbET_RDM34(keys); break;
//Below: RDM3 explicit
    case Coupling::aaaET :
      out = compute_aaaET_RDM34(keys); break;
    case Coupling::bbbET :
      out = compute_bbbET_RDM34(keys); break;
    case Coupling::aabET :
      out = compute_aabET_RDM34(keys); break;
    case Coupling::abbET :
      out = compute_abbET_RDM34(keys); break;
    case Coupling::aETflp :
      out = compute_aETFlip_RDM34(keys); break;
    case Coupling::bETflp :
      out = compute_bETFlip_RDM34(keys); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }
  
  return out;
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aET_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 1: p26B
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // a'b'b'bb
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha}); // a
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B);
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B);
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B);
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0,     0,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactA, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  { //CASE 3', 3'': p27B, 28
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}); // a'b'b
    auto gamma_B1= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aa
    auto gamma_B2= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // b'ab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); 
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B2);
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B2);
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B1);
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    { // E_a'i,bj',ck' 
      int fac = {neleA%2 == 0 ? -1 : 1};
      SMITH::sort_indices<3,2,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
      SMITH::sort_indices<3,2,1,5,0,4, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
      SMITH::sort_indices<3,2,0,4,1,5, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA,     0,     0, nactA,     0, nactA};
      auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
    { // E_a'i',bj',ck 
      int fac = {neleA%2 == 0 ? -1 : 1};
      SMITH::sort_indices<3,4,1,5,0,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); 
      SMITH::sort_indices<3,5,0,4,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); 
      SMITH::sort_indices<3,5,1,4,0,2, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
      SMITH::sort_indices<3,4,0,5,1,2, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA, nactA,     0, nactA,     0,     0};
      auto up  = {nactT, nactT, nactA, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
  }

  { //CASE 5: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha}); //a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'baa
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // b'b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'|a'a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'|a'b'baa 
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // a'|b'b'bba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'a'aaa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'b'baa
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'a'aba
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'b'bba
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA, nactA, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactT, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bET_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 1: p26B
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // b'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // b'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // b'b'b'bb
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta}); // b
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B);
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B);
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B);
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB);
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0,     0,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactA, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  { //CASE 3', 3'': p27B, 28
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}); // b'b'b
    auto gamma_B1= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'ba
    auto gamma_B2= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // b'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); //b'a'a * a'ba
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B2); //b'b'b * b'bb
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B2); //b'a'a * b'bb
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B1); //b'b'b * a'ba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    { // E_a'i,bj',ck' 
      int fac = {neleA%2 == 0 ? -1 : 1};
      SMITH::sort_indices<3,2,1,5,0,4, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a|a'ab
      SMITH::sort_indices<3,2,0,4,1,5, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a|a'ba
      SMITH::sort_indices<3,2,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b|b'bb
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA,     0,     0, nactA,     0, nactA};
      auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
    { // E_a'i',bj',ck 
      int fac = {neleA%2 == 0 ? -1 : 1};
      SMITH::sort_indices<3,5,0,4,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //a'b'a|a'ab
      SMITH::sort_indices<3,4,1,5,0,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //b'b'b|b'bb
      SMITH::sort_indices<3,4,0,5,1,2, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //a'b'a|b'bb
      SMITH::sort_indices<3,5,1,4,0,2, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //b'b'b|a'ab
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA, nactA,     0, nactA,     0,     0};
      auto up  = {nactT, nactT, nactA, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
  }
  { //CASE 5: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta}); //b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'a'aab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'b'bab
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // b'b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'|a'a'aab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'|a'b'bab 
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // b'|b'b'bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'a'aaa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'b'baa
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'a'aba
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'b'bba
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA, nactA, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactT, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaET_RDM34" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 2: p27
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}); // a'a'b'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // aa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'a|aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'a'b'b|aa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj',ck' 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'a|aa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'b|aa
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0, nactA,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  { //CASE 4: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); //a'a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // b'baa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'a'|a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'a'|b'baa 
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck' 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'a'|a'aaa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'a'|b'baa
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA,     0, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bbET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bbET_RDM34" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 2: p27
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // b'b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}); // b'b'b'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'b'a'a|bb
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'b'b|bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj',ck' 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'a|bb
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'b|bb
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0, nactA,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  { //CASE 4: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta}); //b'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // a'abb
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'b'|a'abb
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'b'|b'bbb 
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck' 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'b'|a'abb
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'b'|b'bbb
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA,     0, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abET_RDM34" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 2: p27
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // b'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}); // b'a'b'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // ab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'a'a'a|ab
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'a'b'b|ab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj',ck' 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'a|ab
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'a|ba
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'b|ab
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'b|ba
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0, nactA,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  { //CASE 4: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}); //a'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'aba
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'b'|a'aba
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'b'|b'bba 
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck' 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'b'|a'aba
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a'|a'aab
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'b'|b'bba
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a'|b'bab
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA,     0, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abFlip_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abFlip_RDM34" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 2': p27
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  // b'b'ab 
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}); // a'b
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'b'aa|a'b
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'ab|a'b
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck 
    SMITH::sort_indices<4,2,1,5,0,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'aa|a'b
    SMITH::sort_indices<4,2,1,5,0,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'ab|a'b
    SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'a'bb|a'b (N,M) contribution p41
    SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'ab|a'b (N,M) contribution p41
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA,     0,     0, nactA,     0,     0};
    auto up  = {nactT, nactA, nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  { //CASE 4': p27
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}); // b'a
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'a'ab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); // b'a'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'a|a'a'ab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'a|b'a'bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,b'j',ck' 
    SMITH::sort_indices<3,1,2,4,0,5, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|a'a'ab
    SMITH::sort_indices<3,1,2,4,0,5, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|b'a'bb

    SMITH::sort_indices<5,0,4,3,1,2, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|a'a'ba (N,M) contribution p41
    SMITH::sort_indices<4,0,5,2,1,3, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|a'b'bb (N,M) contribution p41
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA,     0, nactA, nactA,     0, nactA};
    auto up  = {nactT, nactA, nactT, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaaET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaaET_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 3: p26B
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); // a'a'a'
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // aaa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); //a'a'a'|aaa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai',bj',ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0, nactA,     0, nactA,     0, nactA};
    auto up  = {nactA, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bbbET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bbbET_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 3: p26B
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta}); // b'b'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); //b'b'b|bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai',bj',ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB);
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0, nactA,     0, nactA,     0, nactA};
    auto up  = {nactA, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aabET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aabET_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 3: p26B
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha}); // a'b'a
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // aba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); //a'b'a'|aba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai',bj',ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k                                                                                    original order  (cba|ijk)
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'|aba:  (cba|ijk)  
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'|aab:  (bca|ikj) : c->b, k->j
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0, nactA,     0, nactA,     0, nactA};
    auto up  = {nactA, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abbET_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abbET_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 3: p26B
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta}); // b'a'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // bab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // b'a'b'|bab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai',bj',ck' 
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k                                                                                    original order  (cba|ijk)
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'|bab:  (cba|ijk) 
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'|bba:  (bca|ikj) : c->b, k->j
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0, nactA,     0, nactA,     0, nactA};
    auto up  = {nactA, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aETFlip_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aETFlip_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 3': p27B
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}); // a'a'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // b'aa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // a'a'b|b'aa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck' 
    int fac = {neleA%2 == 0 ? -1 : 1};
    SMITH::sort_indices<3,2,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b|b'aa
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA,     0,     0, nactA,     0, nactA};
    auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bETFlip_RDM34(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bETFlip_RDM34" << endl; cout.flush();
  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 3': p27B
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}); // b'b'a
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // a'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // b'b'a|a'bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck' 
    int fac = {neleA%2 == 0 ? -1 : 1};
    SMITH::sort_indices<3,2,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a|a'bb
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA,     0,     0, nactA,     0, nactA};
    auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }

  return make_tuple(out3,out4);
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_diag_RDM34(const array<MonomerKey,4>& keys, const bool subdia) const {
//***************************************************************************************************************
  cout << "DIAG_RDM34" << endl; cout.flush();
  auto& B  = keys[1]; 
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out3 = make_shared<RDM<3>>(nactA+nactB);
  auto out4 = nullptr;

  { //CASE 2'' & 2': p27
    auto gamma_A_aaaa = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'aa
    auto gamma_A_abba = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); //a'b'ba
    auto gamma_A_bbbb = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'bb
    auto gamma_B_aa = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // a'a
    auto gamma_B_bb = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // b'b

    auto rdmA_a  = make_shared<Matrix>(gamma_A_aaaa % gamma_B_aa); //a'a'aa|a'a
    auto rdmA_ab = make_shared<Matrix>(gamma_A_abba % gamma_B_aa); //a'b'ba|a'a
    auto rdmA_b  = make_shared<Matrix>(gamma_A_bbbb % gamma_B_aa); //b'b'bb|a'a
    auto rdmB_a  = make_shared<Matrix>(gamma_A_aaaa % gamma_B_bb); //a'a'aa|b'b
    auto rdmB_ab = make_shared<Matrix>(gamma_A_abba % gamma_B_bb); //a'b'ba|b'b
    auto rdmB_b  = make_shared<Matrix>(gamma_A_bbbb % gamma_B_bb); //b'b'bb|b'b
    
    {
      auto rdmt = rdmA_a->clone();
      // E_a'i',bj,ck
      SMITH::sort_indices<4,5,1,2,0,3, 0,1,  1,1>(rdmA_a->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //a'a'aa|a'a
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|a'a
      SMITH::sort_indices<4,5,0,3,1,2, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmA_b->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //b'b'bb|a'a
      
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmB_a->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //a'a'aa|b'b
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
      SMITH::sort_indices<4,5,0,3,1,2, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|b'b
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmB_b->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //b'b'bb|b'b
      
      if (!subdia) { // <M|op|N> contribution k'j'bc|i'a
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'a'aa|a'a
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|a'a
        SMITH::sort_indices<5,4,3,0,2,1, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'b'bb|a'a

        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmB_a->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'a'aa|b'b
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
        SMITH::sort_indices<5,4,3,0,2,1, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|b'b
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'b'bb|b'b
      }
      
      auto low = {nactA, nactA,     0,     0,     0,     0};
      auto up  = {nactT, nactT, nactA, nactA, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
    }
    {
      auto rdmt = rdmA_a->clone();
      // E_a'i,bj',ck
      SMITH::sort_indices<4,2,1,5,0,3, 0,1, -1,1>(rdmA_a->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //a'a'aa|a'a
      SMITH::sort_indices<4,3,0,5,1,2, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a
      
      SMITH::sort_indices<4,2,1,5,0,3, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
      SMITH::sort_indices<4,2,1,5,0,3, 1,1, -1,1>(rdmB_b->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //b'b'bb|b'b
      
      if (!subdia) { // <M|op|N> contribution k'i'bc|j'a
        SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdmA_a->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //a'a'aa|a'a
        SMITH::sort_indices<5,0,3,4,2,1, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a
        
        SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
        SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdmB_b->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //b'b'bb|b'b
      }
      
      auto low = {nactA,     0,     0, nactA,     0,     0};
      auto up  = {nactT, nactA, nactA, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    }
  }

  { //CASE 4'' & 4': p27
    auto gamma_A_aa = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //a'a
    auto gamma_A_bb = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  //b'b
    auto gamma_B_aaaa = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'aa
    auto gamma_B_abba = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); // a'b'ba
    auto gamma_B_bbbb = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  // b'b'bb

    auto rdmA_a  = make_shared<Matrix>(gamma_A_aa % gamma_B_aaaa); //a'a|a'a'aa
    auto rdmA_ab = make_shared<Matrix>(gamma_A_aa % gamma_B_abba); //a'a|a'b'ba
    auto rdmA_b  = make_shared<Matrix>(gamma_A_aa % gamma_B_bbbb); //a'a|b'b'bb
    auto rdmB_a  = make_shared<Matrix>(gamma_A_bb % gamma_B_aaaa); //b'b|a'a'aa
    auto rdmB_ab = make_shared<Matrix>(gamma_A_bb % gamma_B_abba); //b'b|a'b'ba
    auto rdmB_b  = make_shared<Matrix>(gamma_A_bb % gamma_B_bbbb); //b'b|b'b'bb
    
    {
      auto rdmt = rdmA_a->clone();
      // E_a'i',b'j',ck
      SMITH::sort_indices<3,4,2,5,0,1, 0,1, -1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
      SMITH::sort_indices<3,4,2,5,0,1, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'b'ba
      SMITH::sort_indices<2,5,3,4,0,1, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'a'ab
      SMITH::sort_indices<3,4,2,5,0,1, 1,1, -1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'b'bb

      SMITH::sort_indices<3,4,2,5,0,1, 1,1, -1,1>(rdmB_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'a'aa
      SMITH::sort_indices<3,4,2,5,0,1, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'b'ba
      SMITH::sort_indices<2,5,3,4,0,1, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'a'ab
      SMITH::sort_indices<3,4,2,5,0,1, 1,1, -1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      
      if (!subdia) {
        SMITH::sort_indices<4,3,5,2,1,0, 1,1, -1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
        SMITH::sort_indices<4,3,5,2,1,0, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'b'ba
        SMITH::sort_indices<5,2,4,3,1,0, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'a'ab
        SMITH::sort_indices<4,3,5,2,1,0, 1,1, -1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'b'bb
        
        SMITH::sort_indices<4,3,5,2,1,0, 1,1, -1,1>(rdmB_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'a'aa
        SMITH::sort_indices<4,3,5,2,1,0, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'b'ba
        SMITH::sort_indices<5,2,4,3,1,0, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'a'ab
        SMITH::sort_indices<4,3,5,2,1,0, 1,1, -1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      }
      
      auto low = {nactA, nactA, nactA, nactA,     0,     0};
      auto up  = {nactT, nactT, nactT, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
    }
    {
      auto rdmt = rdmA_a->clone();
      // E_a'i,b'j',ck'
      SMITH::sort_indices<3,1,2,4,0,5, 0,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
      SMITH::sort_indices<2,1,3,4,0,5, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'a'ba
                                                                                                                                   
      SMITH::sort_indices<3,1,2,5,0,4, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'b'ab
      SMITH::sort_indices<3,1,2,4,0,5, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      
      if (!subdia) {
        SMITH::sort_indices<4,0,5,3,1,2, 1,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
        SMITH::sort_indices<5,0,4,3,1,2, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'b'ab
                                                                                                                                     
        SMITH::sort_indices<4,0,5,2,0,3, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'a'ba
        SMITH::sort_indices<4,0,5,3,1,2, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      }
      
      auto low = {nactA, nactA,     0, nactA,     0, nactA};
      auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    }
  }

  return make_tuple(out3,out4);
}



#if 0
//***************************************************************************************************************
double element_RDM3(const int a, const int i, const int b, const int j, const int c, const int k)
//***************************************************************************************************************
  int bar_abc = 0; //monomer abc
  int bar_ijk = 0;
  if(bar_a) ++bar_abc;
  if(bar_b) ++bar_abc;
  if(bar_c) ++bar_abc;
  if(bar_i) ++bar_ijk;
  if(bar_j) ++bar_ijk;
  if(bar_k) ++bar_ijk;
  int bar_tot = bar_abc + bar_ijk;
  switch (bar_tot) {
    case 0 :// monomer A
      break;
    case 1 :// 1 poss.
      break;
    case 2 :// 3 poss.
      break;
    case 3 :// 3 poss.
      break;
    case 4 :// 3 poss.
      break;
    case 5 :// 1 poss.
      break;
    case 6 :// monomer B
      break;
    default :
      throw std::logic_error("3RDM should not reach here.");
      break;
  }
#endif

void ASD_base::symmetrize_RDM34() const {

  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();
  const int nactT = nactA + nactB;  

  //Symmetrize: d(ABAA) note p18B, 19B
  {
    auto low = {0,nactA,0,0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactA); //empty d_ABAA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAA filled
    { //d(AAAB)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactA*nactB); //empty d_AAAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AAAB
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
      SMITH::sort_indices<3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AABA
      auto low = {0,0,nactA,0};
      auto up  = {nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_AABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AABA into d_AABA sector of d
    } 
  }
 
  //Symmetrize: d(ABBB) note p18B, 19B
  {
    auto low = {0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactB); //empty d_ABBB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBB filled
    { //d(BBAB)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactA*nactB); //empty d_BBAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBAB
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
      SMITH::sort_indices<3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBBA
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
    { //d(BAAB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactB); //empty d_BAAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactA); //reorder and fill d_BAAB
      auto low = {nactA,0,0,nactA};
      auto up  = {nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAB into d_BAAB sector of d
    } 
  }

  //Symmetrize: d(AABB) note 19
  { //d(AABB)
    auto low = {0,0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AABB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_AABB
    copy(view.begin(), view.end(), inmat->begin()); //d_AABB filled
    { //d(BBAA)
      auto outmat = make_shared<Matrix>(nactB*nactB*nactA*nactA,1); //empty d_BBAA
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactB, nactB); //reorder and fill d_BBAA
      auto low = {nactA,nactA,0,0};
      auto up  = {nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAA into d_BBAA sector of d
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

}
