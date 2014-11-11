//
// BAGEL - Parallel electron correlation program.
// Filename: asd_rdm4.cc
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

using namespace std;
using namespace bagel;

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>>
ASD_base::couple_blocks_4RDM(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {
//***************************************************************************************************************
  cout << "couple_block_4RDM enetered.." << endl;
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
      out = compute_diag_4RDM(keys, /*subspace diagonal*/false); break;
    case Coupling::aET :
      out = compute_aET_4RDM(keys); break;
    case Coupling::bET :
      out = compute_bET_4RDM(keys); break;
    case Coupling::abFlip :
      out = compute_abFlip_4RDM(keys); break;
    case Coupling::abET :
      out = compute_abET_4RDM(keys); break;
    case Coupling::aaET :
      out = compute_aaET_4RDM(keys); break;
    case Coupling::bbET :
      out = compute_bbET_4RDM(keys); break;
//Below: RDM3 explicit
    case Coupling::aaaET :
      out = compute_aaaET_4RDM(keys); break;
    case Coupling::bbbET :
      out = compute_bbbET_4RDM(keys); break;
    case Coupling::aabET :
      out = compute_aabET_4RDM(keys); break;
    case Coupling::abbET :
      out = compute_abbET_4RDM(keys); break;
    case Coupling::aETflp :
      out = compute_aETFlip_4RDM(keys); break;
    case Coupling::bETflp :
      out = compute_bETFlip_4RDM(keys); break;
// 4RDM
//  case Coupling::a4ET :    
//    out = compute_aaaaET_4RDM(keys); break;
//  case Coupling::b4ET :    
//    out = compute_bbbbET_4RDM(keys); break;
//  case Coupling::a3bET :   
//    out = compute_aaabET_4RDM(keys); break;
//  case Coupling::ab3ET :   
//    out = compute_abbbET_4RDM(keys); break;
//  case Coupling::a2b2ET :  
//    out = compute_aabbET_4RDM(keys); break;
//  case Coupling::doubleFlip :
//    out = compute_doubleFlip_4RDM(keys); break;
//  case Coupling::a2ETflp :
//    out = compute_aaETFlip_4RDM(keys); break;
//  case Coupling::b2ETflp : 
//    out = compute_bbETFlip_4RDM(keys); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }
  
  return out;
}

void ASD_base::initialize_4RDM() {
  cout << "Initialize 4RDM" << endl;
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();

  //E_ai,bj,ck,dl = sum d+c+b+a+ ijkl
  //#of B indices = 0
  fourrdm_.emplace(string("monomerA"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA*nactA*nactA*nactA, 1)); //monomer A
  //# = 1
  fourrdm_.emplace(string("l"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA*nactA*nactA, nactB )); // l
  //# = 2
  fourrdm_.emplace(string("kl"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA*nactA, nactB*nactB )); // kl
  fourrdm_.emplace(string("al"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA*nactA, nactB*nactB )); // al
  fourrdm_.emplace(string("ai"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA*nactA, nactB*nactB )); // ai
  //# = 3
  fourrdm_.emplace(string("jkl"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA, nactB*nactB*nactB )); // jkl
  fourrdm_.emplace(string("akl"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA, nactB*nactB*nactB )); // akl
  fourrdm_.emplace(string("aij"), make_shared<Matrix>(nactA*nactA*nactA*nactA*nactA, nactB*nactB*nactB )); // aij
  //# = 4
  fourrdm_.emplace(string("ijkl"), make_shared<Matrix>(nactA*nactA*nactA*nactA, nactB*nactB*nactB*nactB )); // ijkl
  fourrdm_.emplace(string("ajkl"), make_shared<Matrix>(nactA*nactA*nactA*nactA, nactB*nactB*nactB*nactB )); // ajkl
  fourrdm_.emplace(string("aijk"), make_shared<Matrix>(nactA*nactA*nactA*nactA, nactB*nactB*nactB*nactB )); // aijk
  fourrdm_.emplace(string("bakl"), make_shared<Matrix>(nactA*nactA*nactA*nactA, nactB*nactB*nactB*nactB )); // bakl
  fourrdm_.emplace(string("baij"), make_shared<Matrix>(nactA*nactA*nactA*nactA, nactB*nactB*nactB*nactB )); // baij
  //# = 5
  fourrdm_.emplace(string("aijkl"), make_shared<Matrix>(nactA*nactA*nactA, nactB*nactB*nactB*nactB*nactB )); // aijkl
  fourrdm_.emplace(string("bajkl"), make_shared<Matrix>(nactA*nactA*nactA, nactB*nactB*nactB*nactB*nactB )); // bajkl
  fourrdm_.emplace(string("baijk"), make_shared<Matrix>(nactA*nactA*nactA, nactB*nactB*nactB*nactB*nactB )); // baijk
  //# = 6
  fourrdm_.emplace(string("baijkl"), make_shared<Matrix>(nactA*nactA, nactB*nactB*nactB*nactB*nactB*nactB )); // baijkl
  fourrdm_.emplace(string("cbajkl"), make_shared<Matrix>(nactA*nactA, nactB*nactB*nactB*nactB*nactB*nactB )); // cbajkl
  fourrdm_.emplace(string("cbaijk"), make_shared<Matrix>(nactA*nactA, nactB*nactB*nactB*nactB*nactB*nactB )); // cbaijk
  //# = 7
  fourrdm_.emplace(string("cbaijkl"), make_shared<Matrix>(nactA, nactB*nactB*nactB*nactB*nactB*nactB*nactB )); // cbaijkl
  //# = 8
  fourrdm_.emplace(string("monomerB"), make_shared<Matrix>(nactB*nactB*nactB*nactB*nactB*nactB*nactB*nactB, 1)); // monomer B

  cout << "# of nonredundnat Dimer 4RDM = " << fourrdm_.size() << endl;

}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aET_4RDM" << endl; cout.flush();
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

  //4RDM
  { 
    cout << "4RDM #1" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'a'aaa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta }); //a'b'a'a'aab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta }); //a'b'b'a'abb
    auto gamma_A4 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta }); //a'b'b'b'bbb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha}); // a
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'a'aaa|a
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'b'a'a'aab|a
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // a'b'b'a'abb|a
    auto rdm4 = make_shared<Matrix>(gamma_A4 % gamma_B); // a'b'b'b'bbb|a
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck,dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'a'a'a'aaa|a
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'b'a'a'aab|a
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'a'b'a'aba|a
    SMITH::sort_indices<2,5,1,6,3,4,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'a'a'b'baa|a
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'b'b'a'abb|a
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'b'a'b'bab|a
    SMITH::sort_indices<1,6,3,4,2,5,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'b'b'a'abb|a
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // a'b'b'b'bbb|a
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("l")) += *rdmt;
  }

  { 
    cout << "4RDM #2" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'aa 
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //a'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //a'b'b'bb
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  // b'ab


    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1);  // a'a'a'aa|a'aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1);  // a'b'a'ab|a'aa
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B1);  // a'b'b'bb|a'aa

    auto rdm4 = make_shared<Matrix>(gamma_A1 % gamma_B2);  // a'a'a'aa|b'ab
    auto rdm5 = make_shared<Matrix>(gamma_A2 % gamma_B2);  // a'b'a'ab|b'ab
    auto rdm6 = make_shared<Matrix>(gamma_A3 % gamma_B2);  // a'b'b'bb|b'ab

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck,dl
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<5,6,2,7,1,3,0,4, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a'aa|a'aa
    SMITH::sort_indices<5,6,2,7,1,4,0,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'ba|a'aa
    SMITH::sort_indices<5,6,2,7,0,3,1,4, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'ab|a'aa
    SMITH::sort_indices<5,6,0,7,2,3,1,4, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a'aa|a'aa

    SMITH::sort_indices<5,7,2,6,1,3,1,4, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a'aa|b'ba
    SMITH::sort_indices<5,7,2,6,1,4,0,3, 1,1, -1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'ba|a'aa
    SMITH::sort_indices<5,7,2,6,0,3,1,4, 1,1, -1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'ab|a'aa
    SMITH::sort_indices<5,7,0,6,2,3,1,4, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a'aa|a'aa

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aij")) += *rdmt;
  }
  
  { 
    cout << "4RDM #3" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'aa 
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //a'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //a'b'b'bb

    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  // b'ab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1);  // a'a'a'aa|a'aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1);  // a'b'a'ab|a'aa
    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B2);  // a'b'a'ab|b'ab
    auto rdm4 = make_shared<Matrix>(gamma_A3 % gamma_B2);  // a'b'b'bb|b'ab

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj,ck',dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l
    SMITH::sort_indices<5,3,2,4,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a'aa|a'aa
    SMITH::sort_indices<5,3,1,4,2,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b'ab|a'aa (b<->c)
    SMITH::sort_indices<5,4,2,3,1,7,0,6, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'ba|b'ba (i<->j, k<->l)
    SMITH::sort_indices<5,3,2,4,1,7,0,6, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'bb|b'ba (k<->l)
    SMITH::sort_indices<5,3,2,4,0,6,1,7, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'ba|b'ab (c<->d)
    SMITH::sort_indices<5,3,2,4,0,6,1,7, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'bb|b'ab (c<->d)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("akl")) += *rdmt;
  }

  { 
    cout << "4RDM #4" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  //a'b'b

    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'baa
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); // b'b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'a|a'a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a'a|a'b'baa
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B3); // a'a'a|b'b'bba
 
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B1); // a'b'b|a'a'aaa
    auto rdm5 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'b'b|a'b'baa
    auto rdm6 = make_shared<Matrix>(gamma_A2 % gamma_B3); // a'b'b|b'b'bba

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck',dl
    int fac = {neleA%2 == 0 ? -1 : 1};
    //                  a i b j c k d l
    SMITH::sort_indices<4,5,3,6,1,7,0,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a|a'a'aaa
    SMITH::sort_indices<4,5,3,6,1,7,0,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a|a'b'baa 
    SMITH::sort_indices<3,6,4,5,1,7,0,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a|b'a'aba (a<->b, i<->j)
    SMITH::sort_indices<4,5,3,6,1,7,0,2, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a|b'b'bba

    SMITH::sort_indices<4,5,3,6,0,7,1,2, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b|a'a'aaa (c<->d)
    SMITH::sort_indices<4,5,3,6,0,7,1,2, 1,1, -1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b|a'b'baa (c<->d)
    SMITH::sort_indices<3,6,4,5,0,7,1,2, 1,1, -1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b|b'a'aba (c<->d, a<->b, i<->j)
    SMITH::sort_indices<4,5,3,6,0,7,1,2, 1,1, -1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b|b'b'bba (c<->d)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("baijk")) += *rdmt;
  }

  { 
    cout << "4RDM #5" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  //a'b'b

    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'baa
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); // b'b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'a|a'a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a'a|a'b'baa
 
    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'b'b|a'b'baa
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B3); // a'b'b|b'b'bba

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,b'j',ck',dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a|a'a'aaa
    SMITH::sort_indices<3,2,4,5,1,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a|b'a'baa (a<->b)

    SMITH::sort_indices<4,2,3,6,1,5,0,7, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'b|a'b'aba (j<->k)
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'b|b'b'bba

    SMITH::sort_indices<4,2,3,6,0,7,1,5, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b|a'b'aab (c<->d, l->j, k->l, j->k) (l moved left twice = +1)
    SMITH::sort_indices<4,2,3,5,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b|b'b'bab (c<->d, k<->l)

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("bajkl")) += *rdmt;
  }

  { 
    cout << "4RDM #6" << endl; cout.flush();
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha}); //a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//a'a'a'aaaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//a'a'b'baaa
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//a'b'b'bbaa
    auto gamma_B4 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha});//b'b'b'bbba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'|a'a'a'aaaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'|a'a'b'baaa
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // a'|a'b'b'bbaa
    auto rdm4 = make_shared<Matrix>(gamma_A % gamma_B4); // a'|b'b'b'bbba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',c'k',dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|a'a'a'aaaa
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|a'a'b'baaa
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|a'b'a'abaa (a<->b, i<->j)
    SMITH::sort_indices<1,5,3,6,2,4,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|b'a'a'aaba (a->c, b->a, c->b & i->j, j->k, k->i)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|a'b'b'bbaa
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|b'a'b'baba (b<->c, j<->k)
    SMITH::sort_indices<2,5,1,6,3,4,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|b'b'a'abba (c->a, b->c, a->b &  k->i, i->j, j->k)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // a'|b'b'b'bbba
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("cbaijkl")) += *rdmt;
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bET_4RDM" << endl; cout.flush();
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

  //4RDM 
  { //cf. aET "l"
    cout << "4RDM #1" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'a'a'a'aaa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta }); //b'b'a'a'aab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta }); //b'b'b'a'abb
    auto gamma_A4 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta }); //b'b'b'b'bbb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta}); // b
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'a'a'a'aaa|b
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'a'a'aab|b
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // b'b'b'a'abb|b
    auto rdm4 = make_shared<Matrix>(gamma_A4 % gamma_B); // b'b'b'b'bbb|b
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck,dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'a'a'a'aaa|b
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'b'a'a'aab|b
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'a'b'a'aba|b
    SMITH::sort_indices<2,5,1,6,3,4,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'a'a'b'baa|b
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'b'b'a'abb|b
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'b'a'b'bab|b
    SMITH::sort_indices<1,6,3,4,2,5,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'b'b'a'abb|b
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactA, nactB); // b'b'b'b'bbb|b
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("l")) += *rdmt;
  }

  { //cf aET "aij"
    cout << "4RDM #2" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'a'a'aa 
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //b'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'b'bb
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'ba
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});  // b'bb


    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1);  // b'a'a'aa|a'ba
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1);  // b'b'a'ab|a'ba
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B1);  // b'b'b'bb|a'ba

    auto rdm4 = make_shared<Matrix>(gamma_A1 % gamma_B2);  // b'a'a'aa|b'bb
    auto rdm5 = make_shared<Matrix>(gamma_A2 % gamma_B2);  // b'b'a'ab|b'bb
    auto rdm6 = make_shared<Matrix>(gamma_A3 % gamma_B2);  // b'b'b'bb|b'bb

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck,dl
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l                                                                                                     original order  (dcbkl|aij) => 56271304(=aibjckdl)
    SMITH::sort_indices<5,7,0,6,2,3,1,4, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b'aa|a'ab: -(bdckl|aji) 
    SMITH::sort_indices<5,7,0,6,2,4,1,3, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'ba|a'ab:  (bdclk|aji)
    SMITH::sort_indices<5,7,1,6,2,3,0,4, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'ab|a'ab: -(dbckl|aji)
    SMITH::sort_indices<5,7,2,6,1,3,0,4, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b'bb|a'ab: -(dcbkl|aji)

    SMITH::sort_indices<5,6,1,7,0,3,2,4, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b'aa|b'bb:  (cbdkl|aij)
    SMITH::sort_indices<5,6,0,7,2,4,1,3, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'ba|b'bb:  (bdclk|aij)
    SMITH::sort_indices<5,6,1,7,2,3,0,4, 1,1, -1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'ab|b'bb: -(dbckl|aij)
    SMITH::sort_indices<5,6,2,7,1,3,0,4, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b'bb|b'bb:  (dcbkl|aij)

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aij")) += *rdmt;
  }

  { 
    cout << "4RDM #3" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'a'a'aa 
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //b'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'b'bb

    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'ba
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});  // b'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1);  // b'a'a'aa|a'ba
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1);  // b'b'a'ab|a'ba

    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B2);  // b'b'a'ab|b'bb
    auto rdm4 = make_shared<Matrix>(gamma_A3 % gamma_B2);  // b'b'b'bb|b'bb

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj,ck',dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l                                                                                                     original order  (dcbij|akl)
    SMITH::sort_indices<5,3,2,4,0,6,1,7, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'aa|a'ba: -(cdbij|akl) : c->d
    SMITH::sort_indices<5,3,1,4,0,6,2,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'ab|a'ba:  (cbdij|akl) : d->c->b
    SMITH::sort_indices<5,3,2,4,1,7,0,6, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'aa|a'ab: -(dcbij|alk) : k->l
    SMITH::sort_indices<5,3,1,4,2,7,0,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'ab|a'ab:  (dbcij|alk) : k->l, c->b
    SMITH::sort_indices<5,4,2,3,1,6,0,7, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a'ba|b'bb: -(dcbji|akl) : j->i
    SMITH::sort_indices<5,3,2,4,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b'bb|b'bb:  (dcbij|akl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("akl")) += *rdmt;
  }

  { 
    cout << "4RDM #4" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  //b'b'b

    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'a'aab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'b'bab
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); // b'b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // b'a'a|a'a'aab
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // b'a'a|a'b'bab
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B3); // b'a'a|b'b'bbb
 
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B1); // b'b'b|a'a'aab
    auto rdm5 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'b'b|a'b'bab
    auto rdm6 = make_shared<Matrix>(gamma_A2 % gamma_B3); // b'b'b|b'b'bbb

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck',dl
    int fac = {neleA%2 == 0 ? -1 : 1};
    //                  a i b j c k d l                                                                                                     origianl order  (dcl|baijk)
    SMITH::sort_indices<4,5,3,6,0,7,1,2, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a|a'a'aab: -(cdl|baijk) : d->c
    SMITH::sort_indices<4,5,3,6,0,7,1,2, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a|a'b'bab: -(cdl|baijk) : d->c
    SMITH::sort_indices<3,6,4,5,0,7,1,2, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a|b'a'abb: -(cdl|abjik) : d->c, b->a, i->j
    SMITH::sort_indices<4,5,3,6,1,7,0,2, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a|b'b'bbb: -(cdl|baijk) : d->c

    SMITH::sort_indices<4,5,3,6,1,7,0,2, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b|a'a'aab:  (dcl|baijk)
    SMITH::sort_indices<4,5,3,6,1,7,0,2, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b|a'b'bab:  (dcl|baijk)
    SMITH::sort_indices<3,6,4,5,1,7,0,2, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b|b'a'abb:  (dcl|abjik) : a->b, i->j
    SMITH::sort_indices<4,5,3,6,1,7,0,2, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b|b'b'b'b:  (dcl|baijk)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("baijk")) += *rdmt;
  }

  { 
    cout << "4RDM #5" << endl; cout.flush();
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  //b'b'b

    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'a'aab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'b'bab
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); // b'b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // b'a'a|a'a'aab
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // b'a'a|a'b'bab
 
    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'b'b|a'b'bab
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B3); // b'b'b|b'b'bbb

    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,b'j',ck',dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l                                                                                                     original order  (dci|bajkl)
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'a|a'a'aab   (dci|bajkl)
    SMITH::sort_indices<4,2,3,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a|a'a'aba   (cdi|bajlk) : c->d, k->l
    SMITH::sort_indices<3,2,4,5,1,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'a|b'a'bab  -(dci|abjkl) : a->b
    SMITH::sort_indices<3,2,4,5,1,7,0,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'a|b'a'bba   (dci|abjlk) : a->b, k->l
    SMITH::sort_indices<4,2,3,6,1,5,0,7, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b|a'b'abb  -(dci|abkjl) : j->k
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b|b'b'bbb   (dci|abjkl)

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("bajkl")) += *rdmt;
  }

  { 
    cout << "4RDM #6" << endl; cout.flush();
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta}); //b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});//a'a'a'aaab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});//a'a'b'baab
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});//a'b'b'bbab
    auto gamma_B4 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});//b'b'b'bbbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'|a'a'a'aaab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'|a'a'b'baab
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // b'|a'b'b'bbab
    auto rdm4 = make_shared<Matrix>(gamma_A % gamma_B4); // b'|b'b'b'bbbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',c'k',dl'
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k d l                                                                                                     original order  (d|cbaijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|a'a'a'aaab   (d|cbaijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|a'a'b'baab   (d|cbaijkl)
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|a'b'a'abab   (d|cabjikl) : a->b, i->j
    SMITH::sort_indices<1,6,3,4,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|b'a'a'aabb   (d|acbjkil) : a->b->c, i->j->k
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|a'b'b'bbab   (d|cbaijkl)
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|a'b'b'babb   (d|bcaikjl) : b->c, j->k
    SMITH::sort_indices<1,6,3,4,2,5,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|b'b'a'abbb   (d|acbjkil) : a->b->c, i->j->k
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB, nactB, nactB); // b'|b'b'b'bbbb   (d|cbaijkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("cbaijkl")) += *rdmt;
  }
 
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaET_4RDM" << endl; cout.flush();
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

  //4RDM
  { //p31
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); //a'a'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); //a'a'b'b'bb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // aa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'a'aa|aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'a'b'a'ab|aa
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // a'a'b'b'bb|aa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbaij|kl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'a'aa|aa:  (dcbaij|kl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'a'ab|aa:  (dcbaij|kl)
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'b'ba|aa:  (dcabji|kl) : a->b, i->j
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'b'bb|aa:  (dcbaij|kl)
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("kl")) += *rdmt;
  }

  { //p34B
    cout << "4RDM #2" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});   //a'a'b'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'aaa 
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'baa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'a'a|a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'a'b'b|b'baa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbi|ajkl)
    SMITH::sort_indices<4,3,2,5,1,6,0,7, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'a'a|a'aaa:  (dcbi|ajkl)
    SMITH::sort_indices<4,3,2,5,1,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'b'b|b'baa:  (dcbi|ajkl)
    SMITH::sort_indices<4,3,1,6,2,5,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'a'b|b'aba:  (dbci|akjl) : b->c, j->k
    SMITH::sort_indices<4,3,0,7,2,5,1,6, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'a'b|b'aab:  (bdci|aklj) : b->c->d, j->k->l
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("ajkl")) += *rdmt;
  }

  { //p35
    cout << "4RDM #3" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});   //a'a'b'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'aaa 
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'baa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'a'a|a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a'a'a|b'baa
    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B1); // a'a'b'b|a'aaa
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'a'b'b|b'baa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbl|aijk)
    SMITH::sort_indices<4,5,2,6,1,7,0,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'a'a|a'aaa:  (dcbl|aijk)
    SMITH::sort_indices<4,5,2,6,1,7,0,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'a'a|b'baa:  (dcbl|aijk)
    SMITH::sort_indices<4,5,1,6,0,7,2,3, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'a'b|a'aaa:  (cbdl|aijk) : d->c->b
    SMITH::sort_indices<4,5,1,6,0,7,2,3, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'a'b|b'baa:  (cbdl|aijk) : d->c->b
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aijk")) += *rdmt;
  }
    
  { //p39
    cout << "4RDM #4" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); //a'a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//a'a'aaaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//a'b'baaa
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//b'b'bbaa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'a'|a'a'aaaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'a'|a'b'baaa
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // a'a'|b'b'bbaa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dc|baijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a'|a'a'aaaa:  (dc|baijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a'|a'b'baaa:  (dc|baijkl)
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a'|b'a'abaa:  (dc|abjikl) : b->a, j->i
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a'|b'b'bbaa:  (dc|baijkl)
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("baijkl")) += *rdmt;
  }
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bbET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bbET_4RDM" << endl; cout.flush();
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

  { //p31
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'b'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); //b'b'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); //b'b'b'b'bb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'b'a'a'aa|bb
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'b'a'ab|bb
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // b'b'b'b'bb|bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbaij|kl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'a'aa|bb:  (dcbaij|kl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'a'ab|bb:  (dcbaij|kl)
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'b'ba|bb:  (dcabji|kl) : a->b, i->j
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'b'bb|bb:  (dcbaij|kl)
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("kl")) += *rdmt;
  }

  { //p34B
    cout << "4RDM #2" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //b'b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});   //b'b'b'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //a'abb 
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // b'b'a'a|a'abb
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'b'b'b|b'bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbi|ajkl)
    SMITH::sort_indices<4,3,2,5,1,6,0,7, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'a'a|a'abb:  (dcbi|ajkl)
    SMITH::sort_indices<4,3,2,5,1,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'b'b|b'bbb:  (dcbi|ajkl)
    SMITH::sort_indices<4,3,1,6,2,5,0,7, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'b'a|a'bab:  (dbci|akjl) : b->c, k->j
    SMITH::sort_indices<4,3,1,6,0,7,2,5, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'b'a|a'bba:  (cbdi|aljk) : d->b->c, l->j->k
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("ajkl")) += *rdmt;
  }

  { //p35
    cout << "4RDM #3" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //b'b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});   //b'b'b'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //a'abb 
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // b'b'a'a|a'abb
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // b'b'a'a|b'bbb
    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B1); // b'b'b'b|a'abb
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'b'b'b|b'bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj',ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbl|aijk)
    SMITH::sort_indices<4,5,1,6,0,7,2,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'b'a|a'abb:  (cbdl|aijk) : d->c->b
    SMITH::sort_indices<4,5,1,6,0,7,2,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'b'a|b'bbb:  (cbdl|aijk) : d->c->b
    SMITH::sort_indices<4,5,2,6,1,7,0,3, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'b'b|a'abb:  (dcbl|aijk) 
    SMITH::sort_indices<4,5,2,6,1,7,0,3, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'b'b|b'bbb:  (dcbl|aijk)
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aijk")) += *rdmt;
  }

  { //p39
    cout << "4RDM #4" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta}); //b'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //a'a'aabb
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //a'b'babb
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //b'b'bbbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'b'|a'a'aabb
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'b'|a'b'babb
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // b'b'|b'b'bbbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dc|baijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b'|a'a'aabb:  (dc|baijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b'|a'b'babb:  (dc|baijkl)
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b'|b'a'abbb:  (dc|abjikl) : b->a, j->i
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b'|b'b'bbbb:  (dc|baijkl)
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("baijkl")) += *rdmt;
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abET_4RDM" << endl; cout.flush();
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

  //4RDM
  { //p31
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'b'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //a'b'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //a'b'b'b'bb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // ab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'b'a'a'aa|ab
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'b'b'a'ab|ab
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // a'b'b'b'bb|ab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck',dl'
    //                  a i b j c k d l                                                                                                     original order  (dcbaij|kl)
    SMITH::sort_indices<3,4,2,5,1,7,0,6, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'a'aa|ba:  (dcbaij|lk) : k->l
    SMITH::sort_indices<3,4,2,5,0,6,1,7, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'a'aa|ab:  (cdbaij|kl) : c->d
    SMITH::sort_indices<3,4,2,5,1,7,0,6, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'a'ab|ba:  (dcbaij|lk) : k->l
    SMITH::sort_indices<2,5,3,4,1,7,0,6, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'b'ba|ba:  (dcabji|lk) : k->l, i->j, a->b
    SMITH::sort_indices<3,4,2,5,0,6,1,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'a'ab|ab:  (cdbaij|kl) : c->d
    SMITH::sort_indices<2,5,3,4,0,6,1,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'b'ba|ab:  (cdabji|kl) : c->d, a->b, i->j
    SMITH::sort_indices<3,4,2,5,1,7,0,6, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'b'bb|ba:  (dcbaij|lk) : k->l
    SMITH::sort_indices<3,4,2,5,0,6,1,7, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'b'bb|ba:  (cdbaij|kl) : c->d
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("kl")) += *rdmt;
  }

  { //p34B, 35
    cout << "4RDM #3" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); //b'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}); //b'a'b'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); //a'aba 
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); //b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // b'a'a'a|a'aba
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // b'a'a'a|b'bba
    auto rdm3 = make_shared<Matrix>(gamma_A2 % gamma_B1); // b'a'b'b|a'aba
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'a'b'b|b'bba
    cout << "full gammas" << endl; cout.flush();
    {
      auto rdmt = rdm1->clone();
      // E_a'i,bj',ck',dl' : sign(-1)
      //                  a i b j c k d l                                                                                                     original order  (dcbi|ajkl)
      SMITH::sort_indices<4,3,2,5,1,7,0,6, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'a'a|a'aab  -(dcbi|ajlk) : k->l
      SMITH::sort_indices<4,3,2,5,0,6,1,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'a'a|a'aba  -(cdbi|ajkl) : c->d
      SMITH::sort_indices<4,3,0,6,2,5,1,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'b'a|a'baa  -(bdci|akjl) : b->c->d, j->k
      SMITH::sort_indices<4,3,2,5,1,7,0,6, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'b'b|b'bab  -(dcbi|ajlk) : k->l
      SMITH::sort_indices<4,3,2,5,0,6,1,7, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'b'b|b'bba  -(cdbi|ajkl) : c->d
      SMITH::sort_indices<4,3,1,7,2,5,0,6, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'a'b|b'abb  -(dbci|aklj) : b->c, j->k->l
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("ajkl")) += *rdmt;
    }
    {
      auto rdmt = rdm1->clone();
      // E_a'i',bj',ck',dl : sign(-1)
      //                  a i b j c k d l                                                                                                     original order  (dcbl|aijk)
      SMITH::sort_indices<4,5,2,7,0,6,1,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'a'a|a'aab   (cdbl|aikj) : d->c, j->k
      SMITH::sort_indices<4,5,0,6,1,6,7,3, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'b'a|a'aba   (bcdl|aijk) : b->c->d
      SMITH::sort_indices<4,5,2,7,0,6,1,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'a'a|b'bab   (cdbl|aikj) : d->c, j->k
      SMITH::sort_indices<4,5,0,6,1,6,7,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'b'a|b'bba   (bcdl|aijk) : b->c->d
      SMITH::sort_indices<4,5,1,7,2,6,0,3, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'a'b|a'aab   (dbcl|aikj) : d->c, j->k
      SMITH::sort_indices<4,5,2,6,1,7,0,3, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'b'b|a'aba   (dcbl|aijk) 
      SMITH::sort_indices<4,5,1,7,2,6,0,3, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'a'b|b'bab   (dbcl|aikj) : c->b, j->k
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("aijk")) += *rdmt;
    }

  }

  { //p31
    cout << "4RDM #4" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta} ); //a'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});//a'a'aaba
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});//a'b'baba
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});//b'b'bbba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'b'|a'a'aaba
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'b'|a'b'baba
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // a'b'|b'b'bbba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck',dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dc|baijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'b'|a'a'aaba   (dc|baijkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'a'aaab   (cd|baijlk) : d->c, k->l
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'b'|a'b'baba   (dc|baijkl)
    SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'b'|b'a'abba   (dc|abjikl) : a->b. i->j
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'b'baab   (cd|baijlk) : d->c, k->l
    SMITH::sort_indices<2,5,3,4,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|b'a'abab   (cd|abjilk) : d->c, k->l, a->b. i->j
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'b'|b'b'bbba   (dc|baijkl) 
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|b'b'bbab   (cd|baijlk) : d->c, k->l
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("kl")) += *rdmt;
  }


  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abFlip_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abFlip_4RDM" << endl; cout.flush();
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

  //4RDM
  { //p32
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // b'a'a'aaa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  // b'b'a'aab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  // b'b'b'abb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}); //a'b
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'a'a'aaa|a'b
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'a'aab|a'b
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // b'b'b'abb|a'b
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj,ck,dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dcbijk|al)
    SMITH::sort_indices<6,3,2,4,1,5,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'aaa|a'b   (dcbijk|al)
    SMITH::sort_indices<6,3,2,4,1,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'aab|a'b   (dcbijk|al)
    SMITH::sort_indices<6,3,1,5,2,4,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'aba|a'b   (dbcikj|al) : c->b,  j->k
    SMITH::sort_indices<6,3,2,4,1,5,0,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'abb|a'b   (dcbijk|al)

    //                  a i b j c k d l                                                                                                     original order  (kjibcd|la) : (N,M) conbribution
    SMITH::sort_indices<7,0,3,2,4,1,5,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'aaa|a'b   (ikjbcd|la) : i->j->k
    SMITH::sort_indices<7,1,5,0,3,2,4,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'baa|a'b   (jikcdb|la) : k->j->i, b->c->d
    SMITH::sort_indices<7,1,3,2,5,0,4,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'aba|a'b   (kijbdc|la) : i->j, c->d
    SMITH::sort_indices<7,2,4,1,5,0,3,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bba|a'b   (kjidbc|la) : d->c->b

    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("al")) += *rdmt;
  }

  { //p41B
    cout << "4RDM #2" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});   // b'b'ab
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); //a'a'ab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});  //b'a'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'b'aa|a'a'ab
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'b'ab|b'a'bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,bj,ck,dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dcij|bakl)
    SMITH::sort_indices<5,2,4,3,1,7,0,6, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'aa|a'a'ba  -(dcij|balk) : k->l
    SMITH::sort_indices<5,2,4,3,0,6,1,7, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'aa|a'a'ab  -(cdij|bakl) : d->c
    SMITH::sort_indices<5,2,4,3,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'ab|b'a'bb   (dcij|bakl)
    SMITH::sort_indices<4,3,5,2,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'ba|a'b'bb   (dcji|abkl) : i->j, b->a

    //                  a i b j c k d l                                                                                                     original order  (jicd|lkab) : (N,M) conbribution
    SMITH::sort_indices<7,1,6,0,2,5,3,4, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'aa|a'a'ba  -(jicd|lkba) : a->b
    SMITH::sort_indices<6,0,7,1,2,5,3,4, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'aa|a'a'ab  -(ijcd|lkab) : i->j
    SMITH::sort_indices<6,1,7,0,2,5,3,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'ab|b'a'bb   (jicd|lkab)
    SMITH::sort_indices<6,1,7,0,3,4,2,5, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'ba|a'b'bb   (jidc|klab) : c->d, l->k
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("bakl")) += *rdmt;
  }

  { //p39B
    cout << "4RDM #3" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}); // b'a
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});//a'a'a'aab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});//a'b'a'bab
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});//b'b'a'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'a'|a'a'a'aab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'a'|a'b'a'bab
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // b'a'|b'b'a'bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,b'j',c'k',dl' : sign(-1)
    //                  a i b j c k d l                                                                                                     original order  (di|cbajkl)
    SMITH::sort_indices<4,1,3,5,2,6,0,7, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'a'a'aab  (di|cbajkl)
    SMITH::sort_indices<4,1,3,5,2,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'b'a'bab  (di|cbajkl)
    SMITH::sort_indices<4,1,2,6,3,5,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|b'a'a'abb  (di|bcakjl) : c->b, j->k
    SMITH::sort_indices<4,1,3,5,2,6,0,7, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|b'b'a'bbb  (di|cbajkl)

    //                  a i b j c k d l                                                                                                     original order  (id|lkjabc) : (N,M) conbribution
    SMITH::sort_indices<7,0,5,4,7,3,1,2, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'a'a'baa  (id|lkjbca) : a->b->c
    SMITH::sort_indices<5,0,7,3,6,4,1,2, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'a'b'bba  (id|ljkacb) : k->j, b->c
    SMITH::sort_indices<5,0,6,4,7,3,1,2, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'b'a'bab  (id|lkjabc)
    SMITH::sort_indices<7,0,5,3,7,2,1,4, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a'|a'b'b'bbb  (id|kjlabc) : l->k->j

    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("cbajkl")) += *rdmt;
  }
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaaET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaaET_4RDM" << endl; cout.flush();
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
  
  { //4RDM p32B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // a'a'a'b'b
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //aaa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'a'a|aaa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'a'a'b'b|aaa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_ai,bj',ck',dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a'a'a|aaa   (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a'b'b|aaa   (dcbai|jkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("jkl")) += *rdmt;
  }
  
  { //4RDM p37
    cout << "4RDM #2" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); // a'a'a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aaaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // b'baaa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'a'a'|a'aaaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'a'a'|b'baaa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_ai,bj',ck',dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a'|a'aaaa   (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'a'|b'baaa   (dcb|aijkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aijkl")) += *rdmt;
  }
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bbbET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bbbET_4RDM" << endl; cout.flush();
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
  
  { //4RDM p32B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // b'b'b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // b'b'b'b'b
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'b'b'a'a|bbb
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'b'b'b|bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_ai,bj',ck',dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b'a'a|bbb   (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b'b'b|bbb   (dcbai|jkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("jkl")) += *rdmt;
  }
  
  { //4RDM p37
    cout << "4RDM #2" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta}); // b'b'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // a'abbb
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // b'bbbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'b'b'|a'abbb
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'b'b'|b'bbbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_ai,bj',ck',dl' : sign(+1)
    //                  a i b j c k d l                                                                                                     original order  (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b'|a'abbb   (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'b'|b'bbbb   (dcb|aijkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aijkl")) += *rdmt;
  }
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aabET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aabET_4RDM" << endl; cout.flush();
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

  { //4RDM p32B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // a'b'a'b'b
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); //aba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'b'a'a'a|aba
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'b'a'b'b|aba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_ai,bj',ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'a'a|aba   (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'a'a|aab   (cdbai|jlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b'a'a|baa   (dbcai|kjl) : c-b, j-k

    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a'b'b|aba   (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a'b'b|aab   (cdbai|jlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b'b'b|baa   (dbcai|kjl) : c-b, j-k

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("jkl")) += *rdmt;
  }
  
  { //4RDM p37
    cout << "4RDM #2" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha}); // a'b'a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'aaba
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // b'baba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'b'a'|a'aaba
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'b'a'|b'baba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_a'i',bj',ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a'|a'aaba   (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'a'|a'aaab   (cdb|aijlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'b'|a'abaa   (dbc|aikjl) : c-b, j-k

    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'a'|b'baba   (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'a'|b'baab   (cdb|aijlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'b'|b'bbaa   (dbc|aikjl) : c-b, j-k

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aijkl")) += *rdmt;
  }
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abbET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abbET_4RDM" << endl; cout.flush();
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
  
  { //4RDM p32B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // b'a'b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // b'a'b'b'b
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); //bab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'a'b'a'a|bab
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'a'b'b'b|bab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_ai,bj',ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'a'a|bab   (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'a'a|bba   (cdbai|jlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a'a'a|abb   (dbcai|kjl) : c-b, j-k

    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b'b'b|bab   (dcbai|jkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b'b'b|bba   (cdbai|jlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a'b'b|abb   (dbcai|kjl) : c-b, j-k

    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("jkl")) += *rdmt;
  }
  
  { //4RDM p37
    cout << "4RDM #2" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta}); // b'a'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'abab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // b'bbab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'a'b'|a'abab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'a'b'|b'bbab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_a'i',bj',ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b'|a'abab   (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'b'|a'abba   (cdb|aijlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'a'|a'aabb   (dbc|aikjl) : c-b, j-k

    SMITH::sort_indices<3,4,2,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'a'b'|b'bbab   (dcb|aijkl)
    SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'b'b'|b'bbba   (cdb|aijlk) : c-d, k-l
    SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'a'|b'babb   (dbc|aikjl) : c-b, j-k
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("aijkl")) += *rdmt;
  }
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aETFlip_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aETFlip_4RDM" << endl; cout.flush();
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
  
  //4RDM
  { //p33B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); //a'a'a'ba
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});  //a'a'b'bb
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});//b'aa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'ba|b'aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'a'b'bb|b'aa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_a'i,bj,ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dcbij|akl)
    SMITH::sort_indices<5,3,2,4,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a'ba|b'aa   (dcbij|akl)
    SMITH::sort_indices<5,3,2,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b'bb|b'aa   (dcbij|akl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("akl")) += *rdmt;
  }

  { //p33B
    cout << "4RDM #2" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}); //a'a'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'b'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'b'baa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'a'b|a'b'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'a'b|b'b'baa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_a'i,b'j',ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dci|bajkl)
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'b|a'b'aaa   (dci|bajkl)
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // a'a'b|b'b'baa   (dci|bajkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("bajkl")) += *rdmt;
  }

  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bETFlip_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bETFlip_4RDM" << endl; cout.flush();
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

  //4RDM
  { //p33B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'b'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //b'b'b'ab
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});//a'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'b'a'aa|a'bb
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'b'ab|a'bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_a'i,bj,ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dcbij|akl)
    SMITH::sort_indices<5,3,2,4,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a'aa|a'bb   (dcbij|akl)
    SMITH::sort_indices<5,3,2,4,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b'ab|a'bb   (dcbij|akl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("akl")) += *rdmt;
  }

  { //p33B
    cout << "4RDM #2" << endl;
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}); //b'b'a
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //a'a'abb
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //b'a'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'b'a|a'a'abb
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'b'a|b'a'bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    int fac = {neleA%2 == 0 ? 1 : -1};
    // E_a'i,b'j',ck',dl' 
    //                  a i b j c k d l                                                                                                     original order  (dci|bajkl)
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'a|a'a'abb   (dci|bajkl)
    SMITH::sort_indices<4,2,3,5,1,6,0,7, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB, nactB, nactB); // b'b'a|b'a'bbb   (dci|bajkl)
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    *fourrdm_.at(string("bajkl")) += *rdmt;
  }
  
  return make_tuple(out3,out4);
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_diag_4RDM(const array<MonomerKey,4>& keys, const bool subdia) const {
//***************************************************************************************************************
  cout << "DIAG_4RDM" << endl; cout.flush();
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

  //4RDM
  { //p31B
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'aaa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); //a'b'a'aba
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //b'a'b'bab
    auto gamma_A4 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'b'bbb
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});//a'a
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}); //b'b
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'a'aaa|a'a
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1); // a'b'a'aba|a'a
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B1); // b'a'b'bab|a'a
    auto rdm4 = make_shared<Matrix>(gamma_A4 % gamma_B1); // b'b'b'bbb|a'a

    auto rdm5 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a'a'aaa|b'b
    auto rdm6 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'b'a'aba|b'b
    auto rdm7 = make_shared<Matrix>(gamma_A3 % gamma_B2); // b'a'b'bab|b'b
    auto rdm8 = make_shared<Matrix>(gamma_A4 % gamma_B2); // b'b'b'bbb|b'b
    cout << "full gammas" << endl; cout.flush();

    {
      auto rdmt = rdm1->clone();
      // E_a'i',bj,ck,dl sign(+1)
      //                  a i b j c k d l                                                                                                     original order  (dcbjkl|ai)
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'aaa|a'a   (dcbjkl|ai)
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'aba|a'a   (dcbjkl|ai)
      SMITH::sort_indices<6,7,1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'baa|a'a   (dbckjl|ai) : c->b, j->k
      SMITH::sort_indices<6,7,2,3,0,5,1,4, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'aab|a'a   (cdbjlk|ai) : d->c, k->l
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'bab|a'a   (dcbjkl|ai)
      SMITH::sort_indices<6,7,1,4,2,3,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'abb|a'a   (dbckjl|ai) : c->b, j->k
      SMITH::sort_indices<6,7,2,3,0,5,1,4, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'bba|a'a   (cdbjlk|ai) : d->c, k->l
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bbb|a'a   (dcbjkl|ai)
      
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'aaa|b'b   (dcbjkl|ai)
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'aba|b'b   (dcbjkl|ai)
      SMITH::sort_indices<6,7,1,4,2,3,0,5, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'baa|b'b   (dbckjl|ai) : c->b, j->k
      SMITH::sort_indices<6,7,2,3,0,5,1,4, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'aab|b'b   (cdbjlk|ai) : d->c, k->l
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'bab|b'b   (dcbjkl|ai)
      SMITH::sort_indices<6,7,1,4,2,3,0,5, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'abb|b'b   (dbckjl|ai) : c->b, j->k
      SMITH::sort_indices<6,7,2,3,0,5,1,4, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'bba|b'b   (cdbjlk|ai) : d->c, k->l
      SMITH::sort_indices<6,7,2,3,1,4,0,5, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bbb|b'b   (dcbjkl|ai)
      if (!subdia) { 
        //                  a i b j c k d l                                                                                                     original order  (lkjbcd|ia)
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'aaa|a'a   (lkjbcd|ia)
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'aba|a'a   (lkjbcd|ia)
        SMITH::sort_indices<7,6,4,1,3,2,5,0, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'baa|a'a   (ljkcbd|ia) : k->j, b->c
        SMITH::sort_indices<7,6,3,2,5,0,4,1, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'aab|a'a   (kljbdc|ia) : l->k, c->d
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'bab|a'a   (lkjbcd|ia)
        SMITH::sort_indices<7,6,4,1,3,2,5,0, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'abb|a'a   (ljkcbd|ia) : k->j, b->c
        SMITH::sort_indices<7,6,3,2,5,0,4,1, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'bba|a'a   (lkjbcd|ia) : l->k, c->d
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bbb|a'a   (lkjbcd|ia)

        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'aaa|b'b   (lkjbcd|ia)
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'aba|b'b   (lkjbcd|ia)
        SMITH::sort_indices<7,6,4,1,3,2,5,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'baa|b'b   (ljkcbd|ia) : k->j, b->c
        SMITH::sort_indices<7,6,3,2,5,0,4,1, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'aab|b'b   (kljbdc|ia) : l->k, c->d
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'bab|b'b   (lkjbcd|ia)
        SMITH::sort_indices<7,6,4,1,3,2,5,0, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'abb|b'b   (ljkcbd|ia) : k->j, b->c
        SMITH::sort_indices<7,6,3,2,5,0,4,1, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'bba|b'b   (lkjbcd|ia) : l->k, c->d
        SMITH::sort_indices<7,6,3,2,4,1,5,0, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bbb|b'b   (lkjbcd|ia)
      }
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("ai")) += *rdmt;
    }

    {
      auto rdmt = rdm1->clone();
      // E_a'i,bj,ck,dl' sign(+1)
      //                  a i b j c k d l                                                                                                     original order  (dcbijk|al)
      SMITH::sort_indices<6,3,2,4,1,5,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'aaa|a'a   (dcbijk|al)
      SMITH::sort_indices<6,3,2,5,1,4,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'aab|a'a  -(dcbikj|al) : j->k
      SMITH::sort_indices<6,3,1,4,2,5,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'aba|a'a  -(dbcijk|al) : c->b
      SMITH::sort_indices<6,4,2,3,0,5,1,7, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'abb|a'a   (cdbjik|al) : c->d, j->i

      SMITH::sort_indices<6,4,2,3,0,5,1,7, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'baa|b'b   (dcbijk|al) : c->d, j->i
      SMITH::sort_indices<6,3,1,4,2,5,0,7, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'bab|b'b  -(dbcijk|al) : c->b
      SMITH::sort_indices<6,3,2,5,1,4,0,7, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'bba|b'b   (dcbikj|al) : j->k
      SMITH::sort_indices<6,3,2,4,1,5,0,7, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bbb|b'b   (dcbijk|al)
      if (!subdia) { 
        //                  a i b j c k d l                                                                                                     original order  (kjibcd|la)
        SMITH::sort_indices<7,2,3,1,4,0,5,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'aaa|a'a   (kjibcd|la)
        SMITH::sort_indices<7,2,4,1,3,0,5,6, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'baa|a'a   (kjicbd|la) : c->b
        SMITH::sort_indices<7,2,3,0,4,1,5,6, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'aba|a'a   (jkibcd|la) : j->k
        SMITH::sort_indices<7,1,3,2,5,0,4,6, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'bba|a'a   (kijbdc|la) : j->i, c->d

        SMITH::sort_indices<7,1,3,2,5,0,4,6, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'aab|b'b   (kijbdc|la) : j->i, c->d
        SMITH::sort_indices<7,2,3,0,4,1,5,6, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'bab|b'b   (jkibcd|la) : j->k
        SMITH::sort_indices<7,2,4,1,3,0,5,6, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'abb|b'b   (kjicbd|la) : c->b
        SMITH::sort_indices<7,2,3,1,4,0,5,6, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'bbb|b'b   (kjibcd|la)
      }
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("al")) += *rdmt;
    }
  }

  { //p36, 41
    cout << "4RDM #2" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); //a'b'ba
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'bb
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'aa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); //a'b'ba
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'aa|a'a'aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1); // a'b'ba|a'a'aa
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B1); // b'b'bb|a'a'aa

    auto rdm4 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a'aa|a'b'ba
    auto rdm5 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'b'ba|a'b'ba
    auto rdm6 = make_shared<Matrix>(gamma_A3 % gamma_B2); // b'b'bb|a'b'ba

    auto rdm7 = make_shared<Matrix>(gamma_A1 % gamma_B3); // a'a'aa|b'b'bb
    auto rdm8 = make_shared<Matrix>(gamma_A2 % gamma_B3); // a'b'ba|b'b'bb
    auto rdm9 = make_shared<Matrix>(gamma_A3 % gamma_B3); // b'b'bb|b'b'bb
    cout << "full gammas" << endl; cout.flush();

    {
      auto rdmt = rdm1->clone();
      // E_a'i',b'j',ck,dl sign(+1)
      //                  a i b j c k d l                                                                                                     original order  (dckl|baij)
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|a'a'aa   (dckl|baij)
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|a'a'aa   (dckl|baij)
      SMITH::sort_indices<5,6,4,7,0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|a'a'aa   (cdlk|baij) : c->d, k->l
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|a'a'aa   (dckl|baij)

      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|a'b'ba   (dckl|baij)
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|a'b'ba   (dckl|baij)
      SMITH::sort_indices<5,6,4,7,0,3,1,2, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|a'b'ba   (cdlk|baij) : c->d, k->l
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|a'b'ba   (dckl|baij)

      SMITH::sort_indices<4,7,5,6,1,2,0,3, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|b'a'ab   (dckl|baij) : a->b, i->j
      SMITH::sort_indices<4,7,5,6,1,2,0,3, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|b'a'ab   (dckl|baij) : a->b ,i->j
      SMITH::sort_indices<4,7,5,6,0,3,1,2, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|b'a'ab   (cdlk|baij) : c->d, k->l, a->b, i->j
      SMITH::sort_indices<4,7,5,6,1,2,0,3, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|b'a'ab   (dckl|baij) : a->b, i->j

      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|a'a'aa   (dckl|baij)
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|a'a'aa   (dckl|baij)
      SMITH::sort_indices<5,6,4,7,0,3,1,2, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|a'a'aa   (cdlk|baij) : c->d, k->l
      SMITH::sort_indices<5,6,4,7,1,2,0,3, 1,1,  1,1>(rdm9->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|a'a'aa   (dckl|baij)
      if (!subdia) { 
        //                  a i b j c k d l                                                                                                     original order  (lkcd|jiab)
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|a'a'aa   (lkcd|jiab)
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|a'a'aa   (lkcd|jiab) 
        SMITH::sort_indices<6,5,7,4,3,0,2,1, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|a'a'aa   (kldc|jiab) : k->l, c->d
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|a'a'aa   (lkcd|jiab)

        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|a'b'ba   (lkcd|jiab)
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|a'b'ba   (lkcd|jiab) 
        SMITH::sort_indices<6,5,7,4,3,0,2,1, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|a'b'ba   (kldc|jiab) : k->l, c->d
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|a'b'ba   (lkcd|jiab)

        SMITH::sort_indices<7,4,6,5,2,1,3,0, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|b'a'ab   (lkcd|ijba) : i->j, a->b
        SMITH::sort_indices<7,4,6,5,2,1,3,0, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|b'a'ab   (lkcd|ijba) : i->j, a->b
        SMITH::sort_indices<7,4,6,5,3,0,2,1, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|b'a'ab   (kldc|ijba) : k->l, c->d, i->j, a->b
        SMITH::sort_indices<7,4,6,5,2,1,3,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|b'a'ab   (lkcd|ijba) : i->j, a->b

        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'aa|b'b'bb   (lkcd|jiab)
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'ba|b'b'bb   (lkcd|jiab) 
        SMITH::sort_indices<6,5,7,4,3,0,2,1, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'ab|b'b'bb   (kldc|jiab) : k->l, c->d
        SMITH::sort_indices<6,5,7,4,2,1,3,0, 1,1,  1,1>(rdm9->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'bb|b'b'bb   (lkcd|jiab)
      }
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("baij")) += *rdmt;
    }

    { //p36B, 41B
      auto rdmt = rdm1->clone();
      // E_a'i,b'j,ck',dl' sign(+1)
      //                  a i b j c k d l                                                                                                     original order  (dcij|bakl)
      SMITH::sort_indices<5,2,4,3,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'aa|a'a'aa   (dcij|bakl)
      SMITH::sort_indices<5,2,4,3,1,6,0,7, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'ba|a'b'ba   (dcij|bakl)
      SMITH::sort_indices<4,3,5,2,0,7,1,6, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'ab|b'a'ab   (cdji|ablk) : c->d, i->j, a->b, k->l
      SMITH::sort_indices<4,3,5,2,1,6,0,7, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'ab|b'a'ba   (dcji|abkl) : i->j, b->a
      SMITH::sort_indices<5,2,4,3,0,7,1,6, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'ba|a'b'ab   (cdij|balk) : d->c, k->l
      SMITH::sort_indices<5,2,4,3,1,6,0,7, 1,1,  1,1>(rdm9->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'bb|b'b'bb   (dcij|bakl)
      if (!subdia) { 
        //                  a i b j c k d l                                                                                                     original order  (jicd|lkab)
        SMITH::sort_indices<6,1,7,0,2,5,3,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'a'aa|a'a'aa   (jicd|lkab)
        SMITH::sort_indices<6,1,7,0,2,5,3,4, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'ba|a'b'ba   (jicd|lkab)
        SMITH::sort_indices<7,0,6,1,3,4,2,5, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'ab|b'a'ab   (ijdc|klba) : i->j, c->d, k->l, a->b
        SMITH::sort_indices<6,1,7,0,3,4,2,5, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // a'b'ab|b'a'ba   (jidc|klab) : c->d, k->l
        SMITH::sort_indices<7,0,6,1,2,5,3,4, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'a'ba|a'b'ab   (ijcd|lkba) : i->j, a->b
        SMITH::sort_indices<6,1,7,0,2,5,3,4, 1,1,  1,1>(rdm9->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactA, nactB, nactB); // b'b'bb|b'b'bb   (jicd|lkab)
      }
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("bakl")) += *rdmt;
    }
  }

  { //p39
    cout << "4RDM #1" << endl;
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});//a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}); //b'b
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //a'a'a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); //a'b'a'aba
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  //b'a'b'bab
    auto gamma_B4 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //b'b'b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a|a'a'a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a|a'b'a'aba
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B3); // a'a|b'a'b'bab
    auto rdm4 = make_shared<Matrix>(gamma_A1 % gamma_B4); // a'a|b'b'b'bbb
                                                                 
    auto rdm5 = make_shared<Matrix>(gamma_A2 % gamma_B1); // b'b|a'a'a'aaa
    auto rdm6 = make_shared<Matrix>(gamma_A2 % gamma_B2); // b'b|a'b'a'aba
    auto rdm7 = make_shared<Matrix>(gamma_A2 % gamma_B3); // b'b|b'a'b'bab
    auto rdm8 = make_shared<Matrix>(gamma_A2 % gamma_B4); // b'b|b'b'b'bbb
    cout << "full gammas" << endl; cout.flush();

    {
      auto rdmt = rdm1->clone();
      // E_a'i',bj,ck,dl sign(-1)
      //                  a i b j c k d l                                                                                                     original order  (dl|cbaijk)
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'a'aaa   (dl|cbaijk)
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'a'aba   (dl|cbaijk)
      SMITH::sort_indices<3,6,4,5,2,7,0,1, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'b'baa   (dl|cabjik) : b->a, j->i
      SMITH::sort_indices<4,5,2,7,3,6,0,1, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'a'a'aab   (dl|bcaikj) : b->c, j->k
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'a'b'bab   (dl|cbaijk)
      SMITH::sort_indices<3,6,4,5,2,7,0,1, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'b'a'abb   (dl|cabjik) : b->a, j->i
      SMITH::sort_indices<4,5,2,7,3,6,0,1, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'b'bba   (dl|bcaikj) : b->c, j->k
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'b'b'bbb   (dl|cbaijk)

      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'a'a'aaa   (dl|cbaijk)
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'b'a'aba   (dl|cbaijk)
      SMITH::sort_indices<3,6,4,5,2,7,0,1, 1,1, -1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'a'b'baa   (dl|cabjik) : b->a, j->i
      SMITH::sort_indices<4,5,2,7,3,6,0,1, 1,1, -1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'a'aab   (dl|bcaikj) : b->c, j->k
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'b'bab   (dl|cbaijk)
      SMITH::sort_indices<3,6,4,5,2,7,0,1, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'a'abb   (dl|cabjik) : b->a, j->i
      SMITH::sort_indices<4,5,2,7,3,6,0,1, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'b'b'bba   (dl|bcaikj) : b->c, j->k
      SMITH::sort_indices<4,5,3,6,2,7,0,1, 1,1, -1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'b'bbb   (dl|cbaijk)
      if (!subdia) { 
        //                  a i b j c k d l                                                                                                     original order  (ld|kjiabc)
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'a'aaa   (ld|kjiabc)
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'a'aba   (ld|kjiabc)
        SMITH::sort_indices<6,3,5,4,7,2,1,0, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'b'baa   (ld|kijbac) : j->i, b->a
        SMITH::sort_indices<5,4,7,2,6,3,1,0, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'a'a'aab   (ld|jkiacb) : j->k, b->c
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'a'b'bab   (ld|kjiabc)
        SMITH::sort_indices<6,3,5,4,7,2,1,0, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'b'a'abb   (ld|kijbac) : j->i, b->a
        SMITH::sort_indices<5,4,7,2,6,3,1,0, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'b'bba   (ld|jkiacb) : j->k, b->c
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'b'b'bbb   (ld|kjiabc)

        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm5->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'a'a'aaa   (ld|kjiabc)
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'b'a'aba   (ld|kjiabc)
        SMITH::sort_indices<6,3,5,4,7,2,1,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'a'b'baa   (ld|kijbac) : j->i, b->a
        SMITH::sort_indices<5,4,7,2,6,3,1,0, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'a'aab   (ld|jkiacb) : j->k, b->c
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'b'bab   (ld|kjiabc)
        SMITH::sort_indices<6,3,5,4,7,2,1,0, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'a'abb   (ld|kijbac) : j->i, b->a
        SMITH::sort_indices<5,4,7,2,6,3,1,0, 1,1,  1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'b'b'bba   (ld|jkiacb) : j->k, b->c
        SMITH::sort_indices<5,4,6,3,7,2,1,0, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'b'bbb   (ld|kjiabc)
      }
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("cbaijk")) += *rdmt;
    }

    {
      auto rdmt = rdm1->clone();
      // E_a'i,bj,ck,dl' sign(-1)
      //                  a i b j c k d l                                                                                                     original order  (di|cbajkl)
      SMITH::sort_indices<4,1,3,5,2,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'a'aaa   (di|cbajkl)
      SMITH::sort_indices<4,1,3,6,2,5,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'a'baa  -(di|cbakjl) : k->j
      SMITH::sort_indices<4,1,2,5,3,6,0,7, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'a'a'aba  -(di|bcajkl) : c->b
      SMITH::sort_indices<3,1,4,5,2,7,0,6, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|b'b'a'bba   (di|cabjlk) : a->b, l->k
      SMITH::sort_indices<3,1,4,5,2,7,0,6, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'a|a'a'b'aab   (di|cabjlk) : a->b, l->k
      SMITH::sort_indices<4,1,3,6,2,5,0,7, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'b'abb  -(di|cbakjl) : k->j
      SMITH::sort_indices<4,1,2,5,3,6,0,7, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|a'b'b'bab  -(di|bcajkl) : c->b
      SMITH::sort_indices<4,1,3,5,2,6,0,7, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'b'bbb   (di|cbajkl)
      if (!subdia) { 
        //                  a i b j c k d l                                                                                                     original order  (id|lkjabc)
        SMITH::sort_indices<5,0,6,4,7,3,1,2, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'a'aaa:  (id|lkjabc)
        SMITH::sort_indices<5,0,6,3,7,4,1,2, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'a'b'aba: -(id|ljkabc) : j->k
        SMITH::sort_indices<5,0,7,4,6,3,1,2, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'a'aab: -(id|lkjacb) : c->b
        SMITH::sort_indices<6,0,5,4,7,2,1,3, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // a'a|a'b'b'abb:  (id|kljbac) : l->k, a->b

        SMITH::sort_indices<6,0,5,4,7,2,1,3, 1,1,  1,1>(rdm6->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'a'baa:  (id|kljbac) : l->k, a->b
        SMITH::sort_indices<5,0,7,4,6,3,1,2, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'a'b'bba: -(id|lkjacb) : c->b
        SMITH::sort_indices<5,0,6,3,7,4,1,2, 1,1, -1,1>(rdm7->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'a'bab: -(id|ljkabc) : j->k
        SMITH::sort_indices<5,0,6,4,7,3,1,2, 1,1,  1,1>(rdm8->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB, nactB, nactB); // b'b|b'b'b'bbb:  (id|lkjabc)
      }
      cout << "rearranged" << endl; cout.flush();
      
      *fourrdm_.at(string("cbajkl")) += *rdmt;
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

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaaaET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaaaET_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //4RDM p32B
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); // a'a'a'a'
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //aaaa
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // a'a'a'a'|aaaa
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_ai',bj',ck',dl' : sign(+1)
  //                  a i b j c k d l                                                                                                     original order  (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'a'a'|aaaa   (dcba|ijkl)
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ijkl")) += *rdmt;
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bbbbET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bbbbET_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //4RDM p32B
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta}); // b'b'b'b'
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //bbbb
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // b'b'b'b'|bbbb
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_ai',bj',ck',dl' : sign(+1)
  //                  a i b j c k d l                                                                                                     original order  (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'b'b'|bbbb   (dcba|ijkl)
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ijkl")) += *rdmt;
  
  return make_tuple(out3,out4);
} 

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aabbET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aabbET_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //4RDM p32B
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha}); // a'b'b'a'
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); //abba
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // a'b'b'a'|abba
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_ai',bj',ck',dl' : sign(+1)
  //                  a i b j c k d l                                                                                                     original order  (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'b'a'|abba   (dcba|ijkl)
  SMITH::sort_indices<2,5,3,4,1,6,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'a'b'|baba   (dcab|jikl) : b-a, j-i
  SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'b'a'|abab   (cdba|ijlk) : d-c, k-l
  SMITH::sort_indices<2,5,3,4,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'a'b'|baab   (cdab|jilk) : d-c, b-a, i-j, k-l
  SMITH::sort_indices<2,5,1,6,3,4,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'b'b'|bbaa   (dbac|kijl) : c-b-a, k-j-i
  SMITH::sort_indices<3,4,0,7,2,5,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'a'a'|aabb   (bdca|iklj) : b-c-d, j-k-l
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ijkl")) += *rdmt;
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaabET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaabET_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //4RDM p34
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); //a'b'a'a'
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); //aaba
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // a'b'a'a'|aaba
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_ai',bj',ck',dl' : sign(+1)
  //                  a i b j c k d l                                                                                                     original order  (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'a'a'|aaba   (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'a'a'|aaab   (cdba|ijlk) : c-d, k-l
  SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'b'a'|abaa   (dbca|ikjl) : c-b, j-k
  SMITH::sort_indices<1,6,3,4,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'a'b'|baaa   (dacb|jkil) : a-b-c, i-j-k
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ijkl")) += *rdmt;
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_abbbET_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "abbbET_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();

  //4RDM p34
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::CreateBeta}); //b'a'b'b'
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); //bbab
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // b'a'b'b'|bbab
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_ai',bj',ck',dl' : sign(+1)
  //                  a i b j c k d l                                                                                                     original order  (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'a'b'b'|bbab   (dcba|ijkl)
  SMITH::sort_indices<3,4,2,5,0,7,1,6, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'b'b'b'|bbba   (cdba|ijlk) : c-d, k-l
  SMITH::sort_indices<3,4,1,6,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'a'b'|babb   (dbca|ikjl) : c-b, j-k
  SMITH::sort_indices<1,6,3,4,2,5,0,7, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'b'a'|abbb   (dacb|jkil) : a-b-c, i-j-k
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ijkl")) += *rdmt;
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_aaETFlip_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "aaETFlip_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();
  //4RDM p34B
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}); //a'a'a'b
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'aaa
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // a'a'a'b|b'aaa
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_a'i,bj',ck',dl' : sign(-1)
  //                  a i b j c k d l                                                                                                     original order  (dcbi|ajkl)
  SMITH::sort_indices<4,3,2,5,1,6,0,7, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // a'a'a'b'|baaa   (dcbi|ajkl)
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ajkl")) += *rdmt;
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_bbETFlip_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "bbETFlip_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();
  //4RDM p34B
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}); //b'b'b'a
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); //b'aaa
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // b'b'b'a|a'bbb
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_a'i,bj',ck',dl' : sign(-1)
  //                  a i b j c k d l                                                                                                     original order  (dcbi|ajkl)
  SMITH::sort_indices<4,3,2,5,1,6,0,7, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'b'a'|abbb   (dcbi|ajkl)
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("ajkl")) += *rdmt;
  
  return make_tuple(out3,out4);
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> 
ASD_base::compute_doubleFlip_4RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
  cout << "doubleFlip_4RDM" << endl; cout.flush();
//auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  auto out3 = nullptr;
  auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

//const int neleA = Ap.nelea() + Ap.neleb();
  //4RDM p34B
  cout << "4RDM #1" << endl;
  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); //b'b'aa
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  //a'a'bb
  cout << "partial gammas" << endl; cout.flush();

  auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B); // b'b'aa|a'a'bb
  auto rdmt = rdm1->clone();
  cout << "full gammas" << endl; cout.flush();

  // E_a'i,b'j,ck',dl' : sign(+1)
  //                  a i b j c k d l                                                                                                     original order  (dcij|bakl)
  SMITH::sort_indices<5,2,4,3,1,6,0,7, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'aa|a'a'bb   (dcij|bakl)
  //                  a i b j c k d l                                                                                                     original order  (jicd|lkab) (N,M) contribution
  SMITH::sort_indices<6,1,7,0,2,5,3,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB, nactB, nactB); // b'b'aa|a'a'bb   (jicd|lkab)
  cout << "rearranged" << endl; cout.flush();

  *fourrdm_.at(string("bakl")) += *rdmt;
  
  return make_tuple(out3,out4);
}
