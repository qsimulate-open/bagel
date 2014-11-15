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

double ASD_base::element_3RDM(const int a, const int i, const int b, const int j,const int c, const int k) const {
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  vector<int>  left{c,b,a};
  vector<int>  right{k,j,i};
  vector<bool> bar_left{false,false,false};
  vector<bool> bar_right{false,false,false};
  int nB_left{0};
  //left (creation ops)
  for (int i=0; i<3; i++) {
    int temp = left.at(i) - nactA;
    if(temp >= 0) { //B index
      left.at(i) = temp;
      bar_left.at(i) = true;
      ++nB_left;
    }
  }
  //right (annihilation ops)
  int nB_right{0};
  for (int i=0; i<3; i++) {
    int temp = left.at(i) - nactA;
    if(temp >= 0) { //B index
      right.at(i) = temp;
      bar_right.at(i) = temp;
      ++nB_right;
    }
  }
  //swap left & right
  if(nB_left > nB_right) {
    swap(left,right);
    swap(nB_left, nB_right);
    swap(bar_left, bar_right);
  } 
  //move
  for (int i=1; i<3; i++) { //skip the first element
    if(!bar_right.at(i)) continue; //A index, skip
    //move down B (barred) index
    for (int j=i; j>=1; j--) {
      if(bar_right.at(j-1)) break; //do not interchage with B index
      //swap within right
      swap(right.at(j), right.at(j-1));
      swap(bar_right.at(j), bar_right.at(j-1));
      //swap within left
      swap(left.at(j), left.at(j-1));
      swap(bar_left.at(j), bar_left.at(j-1));
    }
  }
  string bstring;
  //left
  if(bar_left.at(1)) bstring.push_back('c');
  if(bar_left.at(2)) bstring.push_back('b');
  if(bar_left.at(3)) bstring.push_back('a');
  //right
  if(bar_right.at(3)) bstring.push_back('i');
  if(bar_right.at(2)) bstring.push_back('j');
  if(bar_right.at(1)) bstring.push_back('k');

  if(bstring.empty()) bstring = "monomerA";
  if(bstring == "cbaijk") bstring = "monomerB";

  int dim_c = {bar_right.at(1) ? nactB : nactA};
  int dim_b = {bar_right.at(2) ? nactB : nactA};
  int dim_a = {bar_right.at(3) ? nactB : nactA};
  int dim_i = {bar_right.at(3) ? nactB : nactA};
  int dim_j = {bar_right.at(2) ? nactB : nactA};
//int dim_k = {bar_right.at(1) ? nactB : nactA};

  //pointer (a,i,b,j,c,k)
  //      a            i                   b                        j
  int n = left.at(3) + right.at(3)*dim_a + left.at(2)*dim_a*dim_i + right.at(2)*dim_a*dim_i*dim_b +
          left.at(1)*dim_a*dim_i*dim_b*dim_j + right.at(1)*dim_a*dim_i*dim_b*dim_j*dim_c;
  //      c

  //            Matrix             Element
  double* ptr = rdm3_.at(bstring)->element_ptr1(n);
  double out = *ptr;
  delete ptr;

  return out;
}

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
      //                  a i b j c k                                                                                       original order: cbjk|ai
      SMITH::sort_indices<4,5,1,2,0,3, 0,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'a'aa|a'a
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|a'a
      SMITH::sort_indices<4,5,0,3,1,2, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a    bckj|ai
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'b'bb|a'a
      
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmB_a->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //a'a'aa|b'b
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
      SMITH::sort_indices<4,5,0,3,1,2, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|b'b    bckj|ai
      SMITH::sort_indices<4,5,1,2,0,3, 1,1,  1,1>(rdmB_b->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //b'b'bb|b'b
      
      if (!subdia) { // <M|op|N> contribution k'j'bc|i'a
        //                  a i b j c k                                                                                       original order: kjbc|ia
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'a'aa|a'a
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|a'a
        SMITH::sort_indices<5,4,3,0,2,1, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a    jkcb
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'b'bb|a'a

        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmB_a->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'a'aa|b'b
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
        SMITH::sort_indices<5,4,3,0,2,1, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|b'b    jkcb
        SMITH::sort_indices<5,4,2,1,3,0, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'b'bb|b'b
      }
      
      auto low = {nactA, nactA,     0,     0,     0,     0};
      auto up  = {nactT, nactT, nactA, nactA, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
    }

    {
      auto rdmt = rdmA_a->clone();
      // E_a'i,bj',ck (sign: -1)  
      //                  a i b j c k                                                                                       original order : cbik|aj
      SMITH::sort_indices<4,2,1,5,0,3, 0,1, -1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'a'aa|a'a
      SMITH::sort_indices<4,3,0,5,1,2, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a     bcki|aj
      
      SMITH::sort_indices<4,2,1,5,0,3, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
      SMITH::sort_indices<4,2,1,5,0,3, 1,1, -1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'b'bb|b'b
      
      if (!subdia) { // <M|op|N> contribution k'i'bc|j'a
        //                  a i b j c k                                                                                       original order : kibc|ja
        SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdmA_a->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //a'a'aa|a'a
        SMITH::sort_indices<5,0,3,4,2,1, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //b'a'ab|a'a     ikcb
        
        SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); //a'b'ba|b'b
        SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdmB_b->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB);  //b'b'bb|b'b
      }
      //          a     i     b     j     c     k 
      auto low = {nactA,    0,    0,nactA,    0,    0};
      auto up  = {nactT,nactA,nactA,nactT,nactA,nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
    }
  } // Verified

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
      // E_a'i',b'j',ck   (sign: +1)
      //                  a i b j c k                                                                                       original order: ck|baij
      SMITH::sort_indices<3,4,2,5,0,1, 0,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
      SMITH::sort_indices<3,4,2,5,0,1, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'b'ba
      SMITH::sort_indices<2,5,3,4,0,1, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'a'ab    ck|abji
      SMITH::sort_indices<3,4,2,5,0,1, 1,1,  1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'b'bb

      SMITH::sort_indices<3,4,2,5,0,1, 1,1,  1,1>(rdmB_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'a'aa
      SMITH::sort_indices<3,4,2,5,0,1, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'b'ba
      SMITH::sort_indices<2,5,3,4,0,1, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'a'ab    ck|abji
      SMITH::sort_indices<3,4,2,5,0,1, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      
      if (!subdia) {
        //                  a i b j c k                                                                                       original order: kc|jiab
        SMITH::sort_indices<4,3,5,2,1,0, 1,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
        SMITH::sort_indices<4,3,5,2,1,0, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'b'ba
        SMITH::sort_indices<5,2,4,3,1,0, 1,1,  1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'a'ab    kc|ijba
        SMITH::sort_indices<4,3,5,2,1,0, 1,1,  1,1>(rdmA_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'b'bb
        
        SMITH::sort_indices<4,3,5,2,1,0, 1,1,  1,1>(rdmB_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'a'aa  
        SMITH::sort_indices<4,3,5,2,1,0, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'b'ba
        SMITH::sort_indices<5,2,4,3,1,0, 1,1,  1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'a'ab    kc|ijba
        SMITH::sort_indices<4,3,5,2,1,0, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      }
      
      auto low = {nactA, nactA, nactA, nactA,     0,     0};
      auto up  = {nactT, nactT, nactT, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
    }
    {
      auto rdmt = rdmA_a->clone();
      // E_a'i,b'j',ck'   (sign:+1)
      //                  a i b j c k                                                                                       original order: ci|bajk
      SMITH::sort_indices<3,1,2,4,0,5, 0,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
      SMITH::sort_indices<2,1,3,4,0,5, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|b'a'ba   -ci|abjk
                                                                                                                                   
      SMITH::sort_indices<3,1,2,5,0,4, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|a'b'ab   -ci|bakj
      SMITH::sort_indices<3,1,2,4,0,5, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      
      if (!subdia) {
        //                  a i b j c k                                                                                       original order: ic|kjab
        SMITH::sort_indices<4,0,5,3,1,2, 1,1,  1,1>(rdmA_a->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'a'aa
        SMITH::sort_indices<5,0,4,3,1,2, 1,1, -1,1>(rdmA_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //a'a|a'b'ab   -ic|kjba
                                                                                                                                     
        SMITH::sort_indices<4,0,5,2,1,3, 1,1, -1,1>(rdmB_ab->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'a'ba   -ic|jkab
        SMITH::sort_indices<4,0,5,3,1,2, 1,1,  1,1>(rdmB_b->data(),  rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); //b'b|b'b'bb
      }
      //          a     i     b     j     c     k     
      auto low = {nactA,    0,nactA,nactA,    0,nactA};
      auto up  = {nactT,nactA,nactT,nactT,nactA,nactT};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
    }
  } //Verified

  //DEBUG PART
//double sum = 0.0;
  for (int i = 0; i != nactT; ++i)
  for (int j = 0; j != nactT; ++j)
  for (int k = 0; k != nactT; ++k) {
    double elem = out3->element(i,i,j,j,k,k);
  //cout << "CHECK 3RDM(" << i << j << k << ") = " << elem << endl;
    assert(elem > -0.0000000001);
  //sum += elem;
  //if (elem > 1.0e-8) {
  //  cout << "3RDM(" << i << j << k << ") = " << elem << endl;
  //}
  }
//cout << "returning trace.. :" << sum << endl;

  return make_tuple(out3,out4);
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
//auto out4 = nullptr; //make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  //3RDM 
  { //CASE 1: p26B
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  // a'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta});  // a'b'b'bb
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha}); // a
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'aa|a
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'b'a'ab|a
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // a'b'b'bb|a
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck'          even  odd
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k                                                                                    original order: cbaij|k
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); // a'a'a'aa|a
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); // a'b'a'ab|a
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); // a'a'b'ba|a   cabji|k
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); // a'b'b'bb|a
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0,     0,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactA, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified
  
  { //CASE 3', 3'': p27B, 28
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // a'b'b
    auto gamma_B1= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aa
    auto gamma_B2= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});  // b'ab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); // a'a'a| a'aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B1); // a'b'b| a'aa
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B2); // a'a'a| b'ab
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B2); // a'b'b| b'ab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    { // E_a'i,bj',ck'         even odd
      int fac = {neleA%2 == 0 ? -1 : 1};
      //                  a i b j c k                                                                                    original order: cbi|ajk
      SMITH::sort_indices<3,2,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a|a'aa  
      SMITH::sort_indices<3,2,1,5,0,4, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'b'b|b'ba  -cbi|akj
      SMITH::sort_indices<3,2,0,4,1,5, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b|b'ab  -bci|ajk
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA,     0,     0, nactA,     0, nactA};
      auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
    { // E_a'i',bj',ck         even odd
      int fac = {neleA%2 == 0 ? -1 : 1};
      //                  a i b j c k                                                                                    original order: cbk|aij
      SMITH::sort_indices<3,4,1,5,0,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a|a'aa
      SMITH::sort_indices<3,4,0,5,1,2, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b|a'aa  -bck|aij
      SMITH::sort_indices<3,5,1,4,0,2, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'a'a|b'ba  -cbk|aji
      SMITH::sort_indices<3,5,0,4,1,2, 1,1,  1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'b|b'ba   bck|aji
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA, nactA,     0, nactA,     0,     0};
      auto up  = {nactT, nactT, nactA, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
  } //Verified

  { //CASE 5: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha}); //a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'b'baa
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}); // b'b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'|a'a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'|a'b'baa 
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // a'|b'b'bba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck'       even odd
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k                                                                                    original order: c|baijk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'a'aaa 
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'b'baa
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'a'aba   c|abjik
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'b'bba
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA, nactA, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactT, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified

  return make_tuple(out3,nullptr);
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
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});// b'a'a'aa
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // b'b'a'ab
    auto gamma_A3 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); // b'b'b'bb
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta}); // b
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'a'a'aa|b
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'b'a'ab|b
    auto rdm3 = make_shared<Matrix>(gamma_A3 % gamma_B); // b'b'b'bb|b
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj,ck'           even odd
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k                                                                                    original order: cbaij|k 
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); //b'a'a'aa|b
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); //b'b'a'ab|b
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); //b'a'b'ba|b    cabji|k
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactA, nactB); //b'b'b'bb|b
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0,     0,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactA, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified
  
  { //CASE 3', 3'': p27B, 28
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}); // b'a'a
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});  // b'b'b
    auto gamma_B1= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'ba
    auto gamma_B2= gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});   // b'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B1); //b'a'a|a'ba
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B2); //b'b'b|b'bb
    auto rdm3 = make_shared<Matrix>(gamma_A1 % gamma_B2); //b'a'a|b'bb
    auto rdm4 = make_shared<Matrix>(gamma_A2 % gamma_B1); //b'b'b|a'ba
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    { // E_a'i,bj',ck'          even odd
      int fac = {neleA%2 == 0 ? -1 : 1};
      //                  a i b j c k                                                                                       original order: cbi|ajk
      SMITH::sort_indices<3,2,1,5,0,4, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'a'a|a'ab     -     kj
      SMITH::sort_indices<3,2,0,4,1,5, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'b'a|a'ba     -bc
      SMITH::sort_indices<3,2,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'b'b|b'bb
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA,     0,     0, nactA,     0, nactA};
      auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    }
    { // E_a'i',bj',ck         even odd
      int fac = {neleA%2 == 0 ? -1 : 1};
      //                  a i b j c k                                                                                    original order: cbk|aij
      SMITH::sort_indices<3,5,0,4,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //a'b'a|a'ab    bc   ji
      SMITH::sort_indices<3,4,1,5,0,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //b'b'b|b'bb
      SMITH::sort_indices<3,4,0,5,1,2, 1,1, -1,1>(rdm3->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //a'b'a|b'bb   -bc
      SMITH::sort_indices<3,5,1,4,0,2, 1,1, -1,1>(rdm4->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); //b'b'b|a'ab   -     ji
      rdmt->scale(fac);
      cout << "rearranged" << endl; cout.flush();
      
      auto low = {nactA, nactA,     0, nactA,     0,     0};
      auto up  = {nactT, nactT, nactA, nactT, nactA, nactA};
      auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
      copy(rdmt->begin(), rdmt->end(), outv.begin());
      cout << "copied" << endl; cout.flush();
    } // Verified
  }
  { //CASE 5: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta}); //b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'a'aab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'b'bab
    auto gamma_B3 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); // b'b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'|a'a'aab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'|a'b'bab 
    auto rdm3 = make_shared<Matrix>(gamma_A % gamma_B3); // b'|b'b'bbb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',b'j',ck'       even odd
    int fac = {neleA%2 == 0 ? 1 : -1};
    //                  a i b j c k                                                                                       original order: c|baijk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'a'aaa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|a'b'baa
    SMITH::sort_indices<1,4,2,3,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'a'aba        abji
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm3->data(), rdmt->data(), nactA, nactB, nactB, nactB, nactB, nactB); // a'|b'b'bba
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();
    //          a      i      b      j      c      k
    auto low = {nactA, nactA, nactA, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactT, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified

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

    // E_a'i,bj',ck (sign: -1)
    //                  a i b j c k                                                                                       original order: cbik|aj
    SMITH::sort_indices<4,2,1,5,0,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'aa|a'b
    SMITH::sort_indices<4,2,1,5,0,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'ab|a'b
    //(N,M) contribution                                                                                                  original order: kibc|ji
    SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'aa|a'b 
    SMITH::sort_indices<5,1,2,4,3,0, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'ab|a'b 
    cout << "rearranged" << endl; cout.flush();
    //          a      i      b      j      c      k
    auto low = {nactA,     0,     0, nactA,     0,     0};
    auto up  = {nactT, nactA, nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } // Verified
    
  { //CASE 4': p27
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}); // b'a
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // a'a'ab
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}); // b'a'bb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'a|a'a'ab
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'a|b'a'bb
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i,b'j',ck' (sign: +1)
    //                  a i b j c k                                                                                       original order: ci|bajk
    SMITH::sort_indices<3,1,2,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|a'a'ab
    SMITH::sort_indices<3,1,2,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|b'a'bb
    //(N,M) contribution                                                                                                  original order: ic|kjab
    SMITH::sort_indices<5,0,4,3,1,2, 1,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|a'a'ba           ba     
    SMITH::sort_indices<4,0,5,2,1,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a|a'b'bb         jk
    cout << "rearranged" << endl; cout.flush();
    //          a      i      b      j      c      k
    auto low = {nactA,     0, nactA, nactA,     0, nactA};
    auto up  = {nactT, nactA, nactT, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified
  
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
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}); // b'a'b'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}); // ab
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // b'a'a'a|ab
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // b'a'b'b|ab
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj',ck' (sign: +1)
    //                  a i b j c k                                                                                       original order: cbai|jk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'a'a'a|ab
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'a'a|ba      bc  |kj
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'a'b'b|ab
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'b'b'b|ba      bc  |kj
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0, nactA,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified
  
  { //CASE 4: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}); //a'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // a'aba
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}); // b'bba
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'b'|a'aba
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'b'|b'bba 
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck' (sign: +1)
    //                  a i b j c k                                                                                       original order: cb|aijk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'b'|a'aba
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a'|a'aab      bc|  kj
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'b'|b'bba
    SMITH::sort_indices<2,3,0,5,1,4, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'a'|b'bab      bc|  kj
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA,     0, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified

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
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}); // a'a'b'b
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // aa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B); // a'a'a'a|aa
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B); // a'a'b'b|aa
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_ai,bj',ck' (sign: +1)
    //                  a i b j c k                                                                                       original order: cbai|jk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'a'a'a|aa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // a'a'b'b|aa
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0, nactA,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  }//Verified
  
  { //CASE 4: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}); //a'a'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // a'aaa
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}); // b'baa
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // a'a'|a'aaa
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // a'a'|b'baa 
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck' (sign: +1)
    //                  a i b j c k                                                                                       original order: cb|aijk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'a'|a'aaa
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // a'a'|b'baa
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA,     0, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified

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

    // E_ai,bj',ck' (sign: +1)
    //                  a i b j c k                                                                                       original order: cbai|jk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'a'a|bb
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactA, nactB, nactB); // b'b'b'b|bb
    cout << "rearranged" << endl; cout.flush();

    auto low = {    0,     0,     0, nactA,     0, nactA};
    auto up  = {nactA, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified
  
  { //CASE 4: p40B
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta}); //b'b'
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // a'abb
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}); // b'bbb
    cout << "partial gammas" << endl; cout.flush();

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1); // b'b'|a'abb
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2); // b'b'|b'bbb 
    auto rdmt = rdm1->clone();
    cout << "full gammas" << endl; cout.flush();

    // E_a'i',bj',ck' (sign:+1)
    //                  a i b j c k                                                                                       original order: cb|aijk
    SMITH::sort_indices<2,3,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'b'|a'abb
    SMITH::sort_indices<2,3,1,4,0,5, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactB, nactB, nactB, nactB); // b'b'|b'bbb
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA, nactA,     0, nactA,     0, nactA};
    auto up  = {nactT, nactT, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } //Verified

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

    // E_a'i,bj',ck'         even odd
    int fac = {neleA%2 == 0 ? -1 : 1};
    //                  a i b j c k                                                                                       original order: cbi|ajk
    SMITH::sort_indices<3,2,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // a'a'b|b'aa
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA,     0,     0, nactA,     0, nactA};
    auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } // Verified
  
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

    // E_a'i,bj',ck'         even  odd
    int fac = {neleA%2 == 0 ? -1 : 1};
    //                  a i b j c k                                                                                       original order: cbi|ajk
    SMITH::sort_indices<3,2,1,4,0,5, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB, nactB, nactB); // b'b'a|a'bb
    rdmt->scale(fac);
    cout << "rearranged" << endl; cout.flush();

    auto low = {nactA,     0,     0, nactA,     0, nactA};
    auto up  = {nactT, nactA, nactA, nactT, nactA, nactT};
    auto outv = make_rwview(out3->range().slice(low, up), out3->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    cout << "copied" << endl; cout.flush();
  } // Verified

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

  cout << "Symmetrize RDM3.." << endl;
  //#B=1            _            _
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {    0,    0,    0,    0,    0,nactA};
    auto up  = {nactA,nactA,nactA,nactA,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,5),nactB); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(a,a,a,B,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,5),nactB);
      //                  a i c k b j                                           unsorted dimensions
      SMITH::sort_indices<0,1,4,5,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactA, nactA, nactB); 
      auto low = {    0,    0,    0,nactA,    0,    0};
      auto up  = {nactA,nactA,nactA,nactT,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,a,a,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,5),nactB);
      //                  c k a i b j                                           unsorted dimensions
      SMITH::sort_indices<4,5,0,1,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactA, nactA, nactB); 
      auto low = {    0,nactA,    0,    0,    0,    0};
      auto up  = {nactA,nactT,nactA,nactA,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(a,a,a,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,5),nactB);
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactA, nactA, nactB); 
      auto low = {    0,    0,    0,    0,nactA,    0};
      auto up  = {nactA,nactA,nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(a,a,B,a,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,5),nactB);
      //                  i a k c j b                                           unsorted dimensions
      SMITH::sort_indices<1,0,5,4,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactA, nactA, nactB); 
      auto low = {    0,    0,nactA,    0,    0,    0};
      auto up  = {nactA,nactA,nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(B,a,a,a,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,5),nactB);
      //                  k c i a j b                                           unsorted dimensions
      SMITH::sort_indices<5,4,1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactA, nactA, nactB); 
      auto low = {nactA,    0,    0,    0,    0,    0};
      auto up  = {nactT,nactA,nactA,nactA,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#B=2        _   _           __
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {    0,    0,    0,nactA,    0,nactA};
    auto up  = {nactA,nactA,nactA,nactT,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(a,B,a,a,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,    0,    0,    0,nactA};
      auto up  = {nactA,nactT,nactA,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,a,B,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,    0,nactA,    0,    0};
      auto up  = {nactA,nactT,nactA,nactT,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(a,a,B,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactB, nactA, nactB); 
      auto low = {    0,    0,nactA,    0,nactA,    0};
      auto up  = {nactA,nactA,nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(B,a,a,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  j b i a k c                                           unsorted dimensions
      SMITH::sort_indices<3,2,1,0,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,    0,    0,nactA,    0};
      auto up  = {nactT,nactA,nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(B,a,B,a,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  j b k c i a                                           unsorted dimensions
      SMITH::sort_indices<3,2,5,4,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,nactA,    0,    0,    0};
      auto up  = {nactT,nactA,nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#B=2  _     _            _  _  (p43B)
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,    0,    0,nactA,    0,    0};
    auto up  = {nactT,nactA,nactA,nactT,nactA,nactA};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(B,a,a,a,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  a i c k b j                                           unsorted dimensions
      SMITH::sort_indices<0,1,4,5,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactA); 
      auto low = {nactA,    0,    0,    0,    0,nactA};
      auto up  = {nactT,nactA,nactA,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,B,a,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactA); 
      auto low = {    0,nactA,nactA,    0,    0,    0};
      auto up  = {nactA,nactT,nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(a,B,a,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactA); 
      auto low = {    0,nactA,    0,    0,nactA,    0};
      auto up  = {nactA,nactT,nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(a,a,a,B,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  c k b j a i                                           unsorted dimensions
      SMITH::sort_indices<4,5,2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactA); 
      auto low = {    0,    0,    0,nactA,nactA,    0};
      auto up  = {nactA,nactA,nactA,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(a,a,B,a,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  c k a i b j                                           unsorted dimensions
      SMITH::sort_indices<4,5,0,1,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactA); 
      auto low = {    0,    0,nactA,    0,    0,nactA};
      auto up  = {nactA,nactA,nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#=2   _ _                _ _
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,nactA,    0,    0,    0,    0};
    auto up  = {nactT,nactT,nactA,nactA,nactA,nactA};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //1.E(a,a,B,B,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactA, nactA, nactA); 
      auto low = {    0,    0,nactA,nactA,    0,    0};
      auto up  = {nactA,nactA,nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //2.E(a,a,a,a,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactA, nactA, nactA); 
      auto low = {    0,    0,    0,    0,nactA,nactA};
      auto up  = {nactA,nactA,nactA,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#B=3    _   _   _          ___
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {    0,nactA,    0,nactA,    0,nactA};
    auto up  = {nactA,nactT,nactA,nactT,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(B,a,B,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,nactA,    0,nactA,    0};
      auto up  = {nactT,nactA,nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  }

  //#B=3  _ _   _            _ __
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,nactA,    0,nactA,    0,    0};
    auto up  = {nactT,nactT,nactA,nactT,nactA,nactA};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(B,B,a,a,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  a i c k b j                                           unsorted dimensions
      SMITH::sort_indices<0,1,4,5,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {nactA,nactA,    0,    0,    0,nactA};
      auto up  = {nactT,nactT,nactA,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,B,B,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {    0,nactA,nactA,nactA,    0,    0};
      auto up  = {nactA,nactT,nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    }
    { //4.E(a,B,a,a,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {    0,nactA,    0,    0,nactA,nactA};
      auto up  = {nactA,nactT,nactA,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(a,a,a,B,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  c k b j a i                                           unsorted dimensions
      SMITH::sort_indices<4,5,2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {    0,    0,    0,nactA,nactA,nactA};
      auto up  = {nactA,nactA,nactA,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(a,a,B,B,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  c k a i b j                                           unsorted dimensions
      SMITH::sort_indices<4,5,0,1,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {    0,    0,nactA,nactA,    0,nactA};
      auto up  = {nactA,nactA,nactT,nactT,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //7.E(B,B,B,a,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {nactA,nactA,nactA,    0,    0,    0};
      auto up  = {nactT,nactT,nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //8.E(B,B,a,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  i a k c j b                                           unsorted dimensions
      SMITH::sort_indices<1,0,5,4,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {nactA,nactA,    0,    0,nactA,    0};
      auto up  = {nactT,nactT,nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //9.E(B,a,B,B,a,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  j b i a k c                                           unsorted dimensions
      SMITH::sort_indices<3,2,1,0,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {nactA,    0,nactA,nactA,    0,    0};
      auto up  = {nactT,nactA,nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //10.E(B,a,a,a,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  j b k c i a                                           unsorted dimensions
      SMITH::sort_indices<3,2,5,4,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {nactA,    0,    0,    0,nactA,nactA};
      auto up  = {nactT,nactA,nactA,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //11.E(a,a,B,B,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  k c i a j b                                           unsorted dimensions
      SMITH::sort_indices<5,4,1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {    0,    0,nactA,nactA,nactA,    0};
      auto up  = {nactA,nactA,nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //12.E(a,a,B,a,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  k c j b i a                                           unsorted dimensions
      SMITH::sort_indices<5,4,3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactA); 
      auto low = {    0,    0,nactA,    0,nactA,nactA};
      auto up  = {nactA,nactA,nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  }

  //#B=3  _     _   _        _  __ 
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk (p44)
    auto low = {nactA,    0,    0,nactA,    0,nactA};
    auto up  = {nactT,nactA,nactA,nactT,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(a,B,B,a,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,nactA,    0,    0,nactA};
      auto up  = {nactA,nactT,nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,a,B,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,    0,nactA,nactA,    0};
      auto up  = {nactA,nactT,nactA,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(a,B,B,a,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,nactA,    0,nactA,    0};
      auto up  = {nactA,nactT,nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(B,a,a,B,B,a)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  j b i a k c                                           unsorted dimensions
      SMITH::sort_indices<3,2,1,0,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,    0,nactA,nactA,    0};
      auto up  = {nactT,nactA,nactA,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(B,a,B,a,a,B)
      auto outmat = make_shared<Matrix>(pow(nactA,3),pow(nactB,3));
      //                  j b k c i a                                           unsorted dimensions
      SMITH::sort_indices<3,2,5,4,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,nactA,    0,    0,nactA};
      auto up  = {nactT,nactA,nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  }

  //#B=4  _ _   _   _        _ ___
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,nactA,    0,nactA,    0,nactA};
    auto up  = {nactT,nactT,nactA,nactT,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(a,B,B,B,a,B)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,nactA,nactA,    0,nactA};
      auto up  = {nactA,nactT,nactT,nactT,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,a,B,B,B)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactB); 
      auto low = {    0,nactA,    0,nactA,nactA,nactA};
      auto up  = {nactA,nactT,nactA,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(B,B,B,a,B,a)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactB); 
      auto low = {nactA,nactA,nactA,    0,nactA,    0};
      auto up  = {nactT,nactT,nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(B,a,B,B,B,a)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  j b i a k c                                           unsorted dimensions
      SMITH::sort_indices<3,2,1,0,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,nactA,nactA,nactA,    0};
      auto up  = {nactT,nactA,nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(B,a,B,a,B,B)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  j b k c i a                                           unsorted dimensions
      SMITH::sort_indices<3,2,5,4,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactA, nactB, nactA, nactB); 
      auto low = {nactA,    0,nactA,    0,nactA,nactA};
      auto up  = {nactT,nactA,nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#B=4  _   _ _   _      _ _  __
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,    0,nactA,nactA,    0,nactA};
    auto up  = {nactT,nactA,nactT,nactT,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(B,a,a,B,B,B)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  a i c k b j                                           unsorted dimensions
      SMITH::sort_indices<0,1,4,5,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactB, nactB, nactA, nactB); 
      auto low = {nactA,    0,    0,nactA,nactA,nactA};
      auto up  = {nactT,nactA,nactA,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(B,B,B,a,a,B)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  b j a i c k                                           unsorted dimensions
      SMITH::sort_indices<2,3,0,1,4,5, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactB, nactB, nactA, nactB); 
      auto low = {nactA,nactA,nactA,    0,    0,nactA};
      auto up  = {nactT,nactT,nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(B,B,a,B,B,a)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  b j c k a i                                           unsorted dimensions
      SMITH::sort_indices<2,3,4,5,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactB, nactB, nactA, nactB); 
      auto low = {nactA,nactA,    0,nactA,nactA,    0};
      auto up  = {nactT,nactT,nactA,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(a,B,B,B,B,a)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  c k b j a i                                           unsorted dimensions
      SMITH::sort_indices<4,5,2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactB, nactB, nactA, nactB); 
      auto low = {    0,nactA,nactA,nactA,nactA,    0};
      auto up  = {nactA,nactT,nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(a,B,B,a,B,B)
      auto outmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4));
      //                  c k a i b j                                           unsorted dimensions
      SMITH::sort_indices<4,5,0,1,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactA, nactB, nactB, nactA, nactB); 
      auto low = {    0,nactA,nactA,    0,nactA,nactA};
      auto up  = {nactA,nactT,nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#=4   _ _ _ _          _ _ __
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,nactA,nactA,nactA,    0,    0};
    auto up  = {nactT,nactT,nactT,nactT,nactA,nactA};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(nactA*nactA, pow(nactB,4)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //1.E(B,B,a,a,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  a i c k b j                                           unsorted dimensions
      SMITH::sort_indices<0,1,4,5,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactA); 
      auto low = {nactA,nactA,    0,    0,nactA,nactA};
      auto up  = {nactT,nactT,nactA,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //2.E(a,a,B,B,B,B)
      auto outmat = make_shared<Matrix>(pow(nactA,4),nactB*nactB);
      //                  c k a i b j                                           unsorted dimensions
      SMITH::sort_indices<4,5,0,1,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactA); 
      auto low = {    0,    0,nactA,nactA,nactA,nactA};
      auto up  = {nactA,nactA,nactT,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified

  //#B=5  _ _ _ _   _      _ _ ___
  { //1.E(a,i,b,j,c,k) = c'b'a'ijk
    auto low = {nactA,nactA,nactA,nactA,    0,nactA};
    auto up  = {nactT,nactT,nactT,nactT,nactA,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    auto inmat = make_shared<Matrix>(nactA,pow(nactB,5)); //empty
    copy(view.begin(), view.end(), inmat->begin()); //filled
    { //2.E(B,B,a,B,B,B)
      auto outmat = make_shared<Matrix>(nactA,pow(nactB,5));
      //                  a i c k b j                                           unsorted dimensions
      SMITH::sort_indices<0,1,4,5,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactB); 
      auto low = {nactA,nactA,    0,nactA,nactA,nactA};
      auto up  = {nactT,nactT,nactA,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //3.E(a,B,B,B,B,B)
      auto outmat = make_shared<Matrix>(nactA,pow(nactB,5));
      //                  c k a i b j                                           unsorted dimensions
      SMITH::sort_indices<4,5,0,1,2,3, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactB); 
      auto low = {    0,nactA,nactA,nactA,nactA,nactA};
      auto up  = {nactA,nactT,nactT,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //4.E(B,B,B,B,B,a)
      auto outmat = make_shared<Matrix>(nactA,pow(nactB,5));
      //                  i a j b k c                                           unsorted dimensions
      SMITH::sort_indices<1,0,3,2,5,4, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactB); 
      auto low = {nactA,nactA,nactA,nactA,nactA,    0};
      auto up  = {nactT,nactT,nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //5.E(B,B,B,a,B,B)
      auto outmat = make_shared<Matrix>(nactA,pow(nactB,5));
      //                  i a k c j b                                           unsorted dimensions
      SMITH::sort_indices<1,0,5,4,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactB); 
      auto low = {nactA,nactA,nactA,    0,nactA,nactA};
      auto up  = {nactT,nactT,nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
    { //6.E(B,a,B,B,B,B)
      auto outmat = make_shared<Matrix>(nactA,pow(nactB,5));
      //                  k c i a j b                                           unsorted dimensions
      SMITH::sort_indices<5,4,1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactB, nactB, nactB, nactB, nactA, nactB); 
      auto low = {nactA,    0,nactA,nactA,nactA,nactA};
      auto up  = {nactT,nactA,nactT,nactT,nactT,nactT};
      auto outv = btas::make_rwview(threerdm_->range().slice(low,up), threerdm_->storage()); 
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy to threerdm_
    } 
  } //Verified
  cout << "Symmetrize RDM3.. finished" << endl;

}
