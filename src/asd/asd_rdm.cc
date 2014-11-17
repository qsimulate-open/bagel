//
// BAGEL - Parallel electron correlation program.
// Filename: asd_rdm.cc
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
void
ASD_base::debug_RDM() const {
//***************************************************************************************************************
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  const int neleA = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int neleB = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int nelec = neleA+neleB;

  cout << "#of active electrons is A : " << neleA << endl;
  cout << "#of active electrons is B : " << neleB << endl;
  cout << "#of total active electrons: " << nelec << endl;

  auto rdm1A = std::make_shared<RDM<1>>(nactA);
  {
    auto low = {0,0};
    auto up  = {nactA,nactA};
    auto view = btas::make_view(onerdm_->range().slice(low,up), onerdm_->storage());
    copy(view.begin(), view.end(), rdm1A->begin());
  }
  auto rdm1B = std::make_shared<RDM<1>>(nactB);
  {
    auto low = {nactA,nactA};
    auto up  = {nactT,nactT};
    auto view = btas::make_view(onerdm_->range().slice(low,up), onerdm_->storage());
    copy(view.begin(), view.end(), rdm1B->begin());
  }

  auto rdm2A = std::make_shared<RDM<2>>(nactA);
  {
    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage());
    copy(view.begin(), view.end(), rdm2A->begin());
  }
  auto rdm2B = std::make_shared<RDM<2>>(nactB);
  {
    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage());
    copy(view.begin(), view.end(), rdm2B->begin());
  }

//auto rdm3A = std::make_shared<RDM<3>>(nactA);
//{
//  auto low = {0,0,0,0,0,0};
//  auto up  = {nactA,nactA,nactA,nactA,nactA,nactA};
//  auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
//  copy(view.begin(), view.end(), rdm3A->begin());
//}
//auto rdm3B = std::make_shared<RDM<3>>(nactB);
//{
//  auto low = {nactA,nactA,nactA,nactA,nactA,nactA};
//  auto up  = {nactT,nactT,nactT,nactT,nactT,nactT};
//  auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
//  copy(view.begin(), view.end(), rdm3B->begin());
//}

  //1RDM
  { //Monomer A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i) {
      sum += onerdm_->element(i,i);
    }
    std::cout << "1RDM(A)  Trace = " << sum << std::endl;
  }
  { //Monomer B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i) {
      sum += onerdm_->element(i,i);
    }
    std::cout << "1RDM(B)  Trace = " << sum << std::endl;
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i) {
      sum += onerdm_->element(i,i);
    }
    std::cout << "1RDM(AB) Trace = " << sum << std::endl;
  }

  //2RDM check (A) Gamma_ij,kl = <0|E_ij,kl|0> = <0|(k1)'(i2)'(j2)(l1)|0>  1,2 = spin 
  //diagonal: i=j, k=l
  { //Monomer A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j) {
      sum += twordm_->element(i,i,j,j);
    }
    std::cout << "2RDM(A)  Trace = " << sum << std::endl;
  }
  { //Monomer B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j) {
      sum += twordm_->element(i,i,j,j);
    }
    std::cout << "2RDM(B)  Trace = " << sum << std::endl;
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j) {
      sum += twordm_->element(i,i,j,j);
    }
    std::cout << "2RDM(AB) Trace = " << sum << std::endl;
  }

  { //Gamma_ij,kk
    std::cout << "2RDM(A) Partial Trace Sum_k (i,j,k,k)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*rdm1A);
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) {
      debug->element(i,j) -= 1.0/(neleA-1) * twordm_->element(i,j,k,k);
    }
    debug->print(1.0e-3);
  }
  { //Gamma_ij,kk
    std::cout << "2RDM(B) Partial Trace Sum_k (i,j,k,k)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*rdm1B);
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) {
      debug->element(i-nactA,j-nactA) -= 1.0/(neleB-1) * twordm_->element(i,j,k,k);
    }
    debug->print(1.0e-3);
  }
  { //Gamma_ij,kk
    std::cout << "2RDM(AB) Partial Trace Sum_k (i,j,k,k)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*onerdm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) {
      debug->element(i,j) -= 1.0/(nelec-1) * twordm_->element(i,j,k,k);
    }
    debug->print(1.0e-3);
  }

/*
  //construct approx twordm
  cout << "Build approx 2RDM:" << endl;
  for (int i = 0; i != nactA; ++i)
  for (int j = nactA; j != nactT; ++j) {
    approx2rdm_->element(i,i,j,j) = onerdm_->element(i,i) * onerdm_->element(j,j);
    approx2rdm_->element(j,j,i,i) = approx2rdm_->element(i,i,j,j);
    approx2rdm_->element(i,j,j,i) = -0.5*(onerdm_->element(i,i) * onerdm_->element(j,j));// - onerdm_->element(i,i);
    approx2rdm_->element(j,i,i,j) = approx2rdm_->element(i,j,j,i); //-onerdm_->element(i,i) * onerdm_->element(j,j);// - onerdm_->element(j,j);
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j) {
      sum += approx2rdm_->element(i,i,j,j);
    }
    std::cout << "APPROX2RDM(AB) Trace = " << sum << std::endl;
  }
  { //difference
    auto debug = make_shared<RDM<2>>(*twordm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int l = 0; l != nactT; ++l) {
      debug->element(i,j,k,l) -= approx2rdm_->element(i,j,k,l);
    }
    debug->print(1.0e-8);
    auto low = {0,0,0,0};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(debug->range().slice(low,up), debug->storage());
    auto rdm2 = make_shared<Matrix>(nactT*nactT*nactT,nactT,1); 
    copy(view.begin(), view.end(), rdm2->begin());
    cout << "Norm of approx. 2RDM =" << rdm2->norm() << endl;
  }
*/

  //3RDM: Gamma_ij,kl,mn
  std::cout << "-------------- 3RDM --------------" << std::endl;
  { //Trace A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    std::cout << "3RDM Trace (A)  = " << sum << std::endl;
  }
  { //Trace B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    std::cout << "3RDM Trace (B)  = " << sum << std::endl;
  }
  { //Trace AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    std::cout << "3RDM Trace (AB) = " << sum << std::endl;
  }


  { //Gamma_ij,kl,mm : p21
    std::cout << "3RDM Partial Trace Sum_m (i,j,k,l,m,m)" << std::endl;
    auto debug = std::make_shared<RDM<2>>(*twordm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) 
    for (int l = 0; l != nactT; ++l) 
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j,k,l) -= 1.0/(nelec-2) * threerdm_->element(i,j,k,l,m,m);
    }
    debug->print(1.0e-12);
  }
  { //Gamma_ij,mm,kl : p21
    std::cout << "3RDM Partial Trace Sum_m (i,j,m,m,k,l)" << std::endl;
    auto debug = std::make_shared<RDM<2>>(*twordm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) 
    for (int l = 0; l != nactT; ++l) 
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j,k,l) -= 1.0/(nelec-2) * threerdm_->element(i,j,m,m,k,l);
    }
    debug->print(1.0e-12);
  }
  { //Gamma_ij,kl,mm : p21
    std::cout << "3RDM Partial Trace Sum_m (m,m,i,j,k,l)" << std::endl;
    auto debug = std::make_shared<RDM<2>>(*twordm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) 
    for (int l = 0; l != nactT; ++l) 
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j,k,l) -= 1.0/(nelec-2) * threerdm_->element(m,m,i,j,k,l);
    }
    debug->print(1.0e-12);
  }
  { //Gamma_ij,kk,mm : p21
    std::cout << "3RDM Partial Trace Sum_m (i,j,k,k,m,m)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*onerdm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) 
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j) -= 1.0/((nelec-2)*(nelec-1)) * threerdm_->element(i,j,k,k,m,m);
    }
    debug->print(1.0e-12);
  }

  assert(false);
//{ //Gamma_ij,kl,mm : p21
//  std::cout << "3RDM(B) Partial Trace Sum_m (i,j,k,l,m,m)" << std::endl;
//  auto debug = std::make_shared<RDM<2>>(*rdm2B);
//  for (int i = nactA; i != nactT; ++i)
//  for (int j = nactA; j != nactT; ++j)
//  for (int k = nactA; k != nactT; ++k) 
//  for (int l = nactA; l != nactT; ++l) 
//  for (int m = nactA; m != nactT; ++m) {
//    debug->element(i-nactA,j-nactA,k-nactA,l-nactA) -= 1.0/(neleA-2) * threerdm_->element(i,j,k,l,m,m);
//  }
//  debug->print(1.0e-8);
//}
//{ //Gamma_ij,kk,mm : p21
//  std::cout << "3RDM(B) Partial Trace Sum_m (i,j,k,k,m,m)" << std::endl;
//  auto debug = std::make_shared<RDM<1>>(*rdm1B);
//  for (int i = nactA; i != nactT; ++i)
//  for (int j = nactA; j != nactT; ++j)
//  for (int k = nactA; k != nactT; ++k) 
//  for (int m = nactA; m != nactT; ++m) {
//    debug->element(i-nactA,j-nactA) -= 1.0/((neleA-2)*(neleA-1)) * threerdm_->element(i,j,k,k,m,m);
//  }
//  debug->print(1.0e-8);
//}



  assert(false);
}

#if 0

  //4RDM check (A) Gamma_ij,kl,mn,op
  {
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k)
    for (int l = 0; l != nactA; ++l) {
      sum += fourrdm_->element(i,i,j,j,k,k,l,l);
    }
    std::cout << "4RDM Trace = " << sum << std::endl;
  }



  // Checking 4RDM by comparing with 3RDM
  { 
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 1" << std::endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
  { 
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 2" << std::endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
  { 
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 3" << std::endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
  { 
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 4" << std::endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }



  assert(false);

  std::cout << "Monomer A 4RDM print" << std::endl;
  for (int i = 0; i != nactA; ++ i)
  for (int j = 0; j != nactA; ++ j)
  for (int k = 0; k != nactA; ++ k)
  for (int l = 0; l != nactA; ++ l)
  for (int m = 0; m != nactA; ++ m)
  for (int n = 0; n != nactA; ++ n)
  for (int o = 0; o != nactA; ++ o)
  for (int p = 0; p != nactA; ++ p) {
    double elem = fourrdm_->element(i,j,k,l,m,n,o,p);
    elem = std::abs(elem);
    if(elem > 1.0e-8) std::cout << "RDM4(" << i << j << k << l << m << n << o << p << ") " << fourrdm_->element(i,j,k,l,m,n,o,p) << std::endl ;
  }
  assert(false);
#endif

//***************************************************************************************************************
void
ASD_base::debug_energy() const {
//***************************************************************************************************************

  //Energy calculation
  cout << "!@# Energy calculated from RDM:" << endl;
  const int nclosedA = dimer_->active_refs().first->nclosed();
  const int nclosedB = dimer_->active_refs().second->nclosed();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

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

  double  e1 = ddot_(nactT*nactT, int1->element_ptr(0,0), 1, rdm1->element_ptr(0,0), 1);
  cout << "1E energy = " << e1 << endl;

  double e2 = 0.0;
  //AAAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,0>();
    auto int2 = make_shared<Matrix>(nactA*nactA*nactA*nactA,1);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactA); //conver to chemist not.

    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AAAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactA,nactA,1); //empty d_AAAA (note: the dimension specification actually do not matter)
    copy(view.begin(), view.end(), rdm2->begin()); //d_AAAA filled
    cout << "2E energy (AAAA) = " << 0.5 * ddot_(nactA*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //BBBB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,1>();
    auto int2 = make_shared<Matrix>(1,nactB*nactB*nactB*nactB);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactB, nactB, nactB, nactB); //conver to chemist not.

    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBBB sector of d
    auto rdm2 = make_shared<Matrix>(1,nactB*nactB*nactB*nactB); //empty d_BBBB
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBBB filled
    cout << "2E energy (BBBB) = " << 0.5 * ddot_(nactB*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactB*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  cout << "(3A,1B) part:" << endl;
  //AAAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,1>(); // <pq|rs'> in (pqr,s') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactA*nactB,1);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.
    auto low = {    0,    0,    0,nactA};
    auto up  = {nactA,nactA,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AAAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactA*nactB,1); //empty d_AAAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_AAAB filled
    cout << "2E energy (AAAB) = " << 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //AABA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,0>(); // <pq'|rs> in (prs,q') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactA,1);
    SMITH::sort_indices<0,1,3,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [ij|k'l]

    auto low = {    0,    0,nactA,    0};
    auto up  = {nactA,nactA,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AABA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactA,1); //empty d_AABA
    copy(view.begin(), view.end(), rdm2->begin()); //d_AABA filled
    cout << "2E energy (AABA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //BAAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,0>(); // <p'q|rs> in (qrs,p') format 
    auto int2 = make_shared<Matrix>(nactB,nactA*nactA*nactA);
    SMITH::sort_indices<3,1,0,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.

    auto low = {nactA,    0,    0,    0};
    auto up  = {nactT,nactA,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAA sector of d
    auto rdm2 = make_shared<Matrix>(nactB,nactA*nactA*nactA); //empty d_BAAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BAAA filled
    cout << "2E energy (BAAA) = " << 0.5 * ddot_(nactB*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,0>(); // <pq|r's> in (pqs,r') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactA*nactA,1);
    SMITH::sort_indices<0,3,1,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [pr'|qs]=[ij'|kl]

    auto low = {    0,nactA,    0,    0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactA*nactA,1); //empty d_ABAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABAA filled
    cout << "2E energy (ABAA) = " << 0.5 * ddot_(nactA*nactB*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }



  cout << "(1A,3B) part:" << endl;
  //ABBB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,1>(); // <pq'|r's'> in (p,q'r's') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [pr'|q's']=[ij'|k'l']

    auto low = {    0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_ABBB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (ABBB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BABB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,1>(); // <p'q'|rs'> in (r,p'q's') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<1,0,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r|q's']=[i'j|k'l']

    auto low = {nactA,    0,nactA,nactA};
    auto up  = {nactT,nactA,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BABB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BABB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (BABB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,1>(); // <p'q|r's'> in (q,p'r's') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<1,2,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|qs']=[i'j'|kl']

    auto low = {nactA,nactA,    0,nactA};
    auto up  = {nactT,nactT,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BBAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (BBAB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBBA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,0>(); // <p'q'|r's> in (s,p'q'r') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<1,3,2,0, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|q's]=[i'j'|k'l]

    auto low = {nactA,nactA,nactA,    0};
    auto up  = {nactT,nactT,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBBA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BBBA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBBA filled
    cout << "2E energy (BBBA) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  cout << "(2A,2B) part:" << endl;
  //AABB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,1>(); // <pq'|rs'> in (pr,q's') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<0,1,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr|q's']=[ij|k'l']

    auto low = {    0,    0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AABB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_AABB
    copy(view.begin(), view.end(), rdm2->begin()); //d_AABB filled
    cout << "2E energy (AABB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,0>(); // <p'q|r's> in (qs,p'r') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r'|qs]=[i'j'|kl]

    auto low = {nactA,nactA,    0,    0};
    auto up  = {nactT,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BBAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBAA filled
    cout << "2E energy (BBAA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,1>(); // <pq|r's'> in (pq,r's') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|qs']=[ij'|kl']

    auto low = {    0,nactA,    0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_ABAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABAB filled
    cout << "2E energy (ABAB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BABA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,0>(); // <p'q'|rs> in (rs,p'q') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<2,0,3,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|q's]=[i'j|k'l]

    auto low = {nactA,    0,nactA,    0};
    auto up  = {nactT,nactA,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BABA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BABA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BABA filled
    cout << "2E energy (BABA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABBA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,0>(); // <pq'|r's> in (ps,q'r') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<0,3,2,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|q's]=[ij'|k'l]

    auto low = {    0,nactA,nactA,    0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_ABBA
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBA filled
    cout << "2E energy (ABBA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BAAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,1>(); // <p'q|rs'> in (qr,p's') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<2,1,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|qs']=[i'j|kl']

    auto low = {nactA,    0,    0,nactA};
    auto up  = {nactT,nactA,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BAAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_BAAB filled
    cout << "2E energy (BAAB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  
  //Energy print
  cout << "nuclear repulsion= " << dimer_->sref()->geom()->nuclear_repulsion() << endl;
  cout << "core energy      = " << jop_->core_energy() << endl;
  cout << "nuc + core energy= " << dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy() << endl;
  cout << "1E energy = " <<  e1 << endl;
  cout << "2E energy = " <<  e2 << endl;
  cout << "Total energy = " << dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy() + e1 + e2 << endl;
}

