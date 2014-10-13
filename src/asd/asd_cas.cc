//
// BAGEL - Parallel electron correlation program.
// Filename: asd_cas.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/asd_cas.h>
#include <src/fci/hztasks.h>
#include <src/fci/prop1etask.h>

#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

shared_ptr<Dvec> ASD_CAS::form_sigma(shared_ptr<const Dvec> ccvec, std::shared_ptr<const MOFile> jop) const {
  const int nstates = ccvec->ij();

  shared_ptr<const Determinants> det = ccvec->det();
  auto int_det = det->remalpha()->rembeta();

  auto sigmavec = make_shared<Dvec>(det, nstates);
  sigmavec->zero();

  const int norb = det->norb();

  shared_ptr<Matrix> h1 = jop->mo1e()->matrix();

  auto h2 = make_shared<Matrix>(norb*norb, norb*norb);
  double* h2_ptr = h2->data();
  for (int i = 0, ijkl = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      for (int k = 0; k < norb; ++k) {
        for (int l = 0; l < norb; ++l, ++ijkl) {
          h2_ptr[ijkl] = jop->mo2e_hz(l,k,j,i) - jop->mo2e_hz(k,l,j,i);
        }
      }
    }
  }

  const int ij = norb * norb;
  auto d = make_shared<Dvec>(int_det, ij);
  auto e = make_shared<Dvec>(int_det, ij);

  for (int istate = 0; istate < nstates; ++istate) {
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    sigma_aa(cc, sigma, h1->data(), h2->data());

    auto cc_trans = cc->transpose();
    auto sg_trans = make_shared<Civec>(cc_trans->det());

    // sigma_bb
    sigma_aa(cc_trans, sg_trans, h1->data(), h2->data());

    sigma->ax_plus_y(1.0, *sg_trans->transpose());

    d->zero();

    sigma_2ab_1(cc, d);
    sigma_2ab_2(d, e, jop->mo2e_ptr());
    sigma_2ab_3(sigma, e);
  }

  return sigmavec;
}

shared_ptr<Dvec> ASD_CAS::form_sigma_1e(shared_ptr<const Dvec> ccvec, const double* modata) const {
  const int nstate = ccvec->ij();
  shared_ptr<const Determinants> det = ccvec->det();

  const int lbs = det->lenb();
  const int las = det->lena();

  shared_ptr<const Dvec> cc_trans = ccvec->spinflip();
  shared_ptr<const Determinants> det_trans = cc_trans->det();

  auto sigma = make_shared<Dvec>(det, nstate);
  auto sg_trans = make_shared<Dvec>(det_trans, nstate);

  TaskQueue<Prop1eTask> tasks((det->lena() + det_trans->lenb()) * nstate);

  for (int istate = 0; istate < nstate; ++istate) {
    double* target = sigma->data(istate)->data();
    for (auto& a : det->string_bits_a()) {
      tasks.emplace_back(ccvec->data(istate), a, target, modata);
      target += lbs;
    }

    target = sg_trans->data(istate)->data();
    for (auto& a : det_trans->string_bits_a()) {
      tasks.emplace_back(cc_trans->data(istate), a, target, modata);
      target += las;
    }
  }

  tasks.compute();

  sigma->ax_plus_y(1.0, *sg_trans->spinflip());

  return sigma;
}


void ASD_CAS::sigma_aa(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, const double* const h1, const double* const h2) const {
  assert(*cc->det() == *sigma->det());

  shared_ptr<const Determinants> det = cc->det();
  const int lb = cc->lenb();

  TaskQueue<HZTaskAA<double>> tasks(det->lena());

  double* target = sigma->data();
  for (auto& a : det->string_bits_a()) {
    tasks.emplace_back(cc, a, target, h1, h2);
    target += lb;
  }

  tasks.compute();
}


void ASD_CAS::sigma_2ab_1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {

  shared_ptr<const Determinants> base_det = cc->det();
  shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta();

  const int norb = base_det->norb();
  const int lbs = base_det->lenb();
  const double* source_base = cc->data();

  TaskQueue<HZTaskAB1<double>> tasks(norb*norb);

  for (int k = 0; k < norb; ++k) {
    for (int l = 0; l < norb; ++l) {
      double* target_base = d->data(k*norb + l)->data();
      tasks.emplace_back(int_det, lbs, source_base, target_base, k, l);
    }
  }

  tasks.compute();
}

void ASD_CAS::sigma_2ab_2(shared_ptr<Dvec> d, shared_ptr<Dvec> e, const double* mo2e_ptr) const {
  const int lenab = d->lena() * d->lenb();
  const int ij = d->ij();

  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, mo2e_ptr, ij, 0.0, e->data(), lenab);
}

void ASD_CAS::sigma_2ab_3(shared_ptr<Civec> sigma, shared_ptr<Dvec> e) const {
  shared_ptr<const Determinants> base_det = sigma->det();
  shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta();

  const int norb = base_det->norb();
  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  double* target_base = sigma->data();

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      const double* source_base = e->data(i*norb + j)->data();
      for (auto& aiter : int_det->phiupa(i)) {
        double *target = target_base + aiter.target*lbt;
        const double *source = source_base + aiter.source*lbs;
        for (auto& biter : int_det->phiupb(j)) {
          const double sign = static_cast<double>(aiter.sign * biter.sign);
          target[biter.target] += sign * source[biter.source];
        }
      }
    }
  }
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>
ASD_CAS::compute_rdm34_monomer (pair<int,int> offset, array<Dvec,4>& fourvecs) const {
// taken from fci/fci_rdm.cc
//***************************************************************************************************************

  std::cout << "ASD_CAS: compute_rdm34_monomer called: " << get<0>(offset) << std::endl;
  
  //Dvec decomposition
  auto& A  = fourvecs[0]; // <I'>
  auto& Ap = fourvecs[1]; // |I>

  auto& B  = fourvecs[2]; // <J'|
  auto& Bp = fourvecs[3]; // |J>

  //Offsets
  const int ioff = std::get<0>(offset);
  const int joff = std::get<1>(offset);
  assert(ioff == joff); //diagonal only at the moment
  
  const int nactA = dimer_->active_refs().first->nact(); //dimer_ (ASD_base)
  const int nactB = dimer_->active_refs().second->nact();
  const int nactT = nactA+nactB;

  auto rdm3A = std::make_shared<RDM<3>>(nactA);
//auto rdm4A = std::make_shared<RDM<4>>(nactA);
  auto rdm3B = std::make_shared<RDM<3>>(nactB);
//auto rdm4B = std::make_shared<RDM<4>>(nactB);

  const int nstA  = A.ij();
  const int nstAp = Ap.ij();
  std::cout << "<I'| x |I> = " << nstA << " x " << nstAp << std::endl;

  const int nstB  = B.ij();
  const int nstBp = Bp.ij();
  std::cout << "<J'| x |J> = " << nstA << " x " << nstAp << std::endl;

  assert(nstA == nstB && nstAp == nstBp);

  //currently, Diagonal subspace only
  assert(nstA == nstAp);

  //ground state only
  //const int istate = 0;

  //MonomerA
  for (int i = 0; i != nstA; ++i) {//<I'|
    for (int ip = 0; ip != nstAp; ++ip) {// |I>
  
      shared_ptr<RDM<3>> r3;
      shared_ptr<RDM<4>> r4;
      tie(r3,r4) = compute_rdm34_from_civec(A.data(i), Ap.data(ip)); // <I'|E(op)|I>

      double csum = 0.0; //coeff sum
      for (int j = 0; j != nstB; ++j) { // delta_J'J
        const int ij  = i  + (j*nstA);
        const int ijp = ip + (j*nstAp);
        csum += adiabats_->element(ioff+ij,0) * adiabats_->element(joff+ijp,0);
      } //I'I

      r3->scale(csum);
      *rdm3A += *r3;
//    r4->scale(csum);
//    *rdm4A += *r4;
    } //|I>
  } //<I'|
  //MonomerB
  for (int j = 0; j != nstB; ++j) {//<J'|
    for (int jp = 0; jp != nstBp; ++jp) {// |J>
  
      shared_ptr<RDM<3>> r3;
      shared_ptr<RDM<4>> r4;
      tie(r3,r4) = compute_rdm34_from_civec(B.data(j), Bp.data(jp)); // <J'|E(op)|J>

      double csum = 0.0; //coeff sum
      for (int i = 0; i != nstA; ++i) { // delta_I'I
        const int ij  = i + (j*nstA);
        const int ijp = i + (jp*nstAp);
        csum += adiabats_->element(ioff+ij,0) * adiabats_->element(joff+ijp,0);
      } //I'I

      r3->scale(csum);
      *rdm3B += *r3;
//    r4->scale(csum);
//    *rdm4B += *r4;
    } //|J>
  } //<J'|
  //END NEW

  cout << "partial rdm complete" << endl; cout.flush();
  auto out3 = std::make_shared<RDM<3>>(nactA+nactB);
  out3->zero();
  {
    //Monomer A
    auto low = {    0,    0,    0,    0,    0,    0};
    auto up  = {nactA,nactA,nactA,nactA,nactA,nactA};
    auto outv = make_rwview(out3->range().slice(low,up), out3->storage());
    copy(rdm3A->begin(), rdm3A->end(), outv.begin());
  }
  {
    //Monomer B
    auto low = {nactA,nactA,nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT,nactT,nactT};
    auto outv = make_rwview(out3->range().slice(low,up), out3->storage());
    copy(rdm3B->begin(), rdm3B->end(), outv.begin());
  }
  cout << "rdm3 copy complete" << endl; cout.flush();
  shared_ptr<RDM<4>> out4 = nullptr;
//auto out4 = std::make_shared<RDM<4>>(nactA+nactB);
//out4->zero();
//{
//  //Monomer A
//  auto low = {    0,    0,    0,    0,    0,    0,    0,    0};
//  auto up  = {nactA,nactA,nactA,nactA,nactA,nactA,nactA,nactA};
//  auto outv = make_rwview(out4->range().slice(low,up), out4->storage());
//  copy(rdm4A->begin(), rdm4A->end(), outv.begin());
//}
//{
//  //Monomer B
//  auto low = {nactA,nactA,nactA,nactA,nactA,nactA,nactA,nactA};
//  auto up  = {nactT,nactT,nactT,nactT,nactT,nactT,nactT,nactT};
//  auto outv = make_rwview(out4->range().slice(low,up), out4->storage());
//  copy(rdm4B->begin(), rdm4B->end(), outv.begin());
//}
//cout << "rdm4 copy complete" << endl; cout.flush();

  return make_tuple(out3, out4);

}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>>
ASD_CAS::compute_rdm34_from_civec (shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {
//***************************************************************************************************************
  cout << "ASD_CAS::compute_rdm34_from_civec" << endl; cout.flush();
  //ADDED
  const int norb = cbra->det()->norb();
  //END

  auto rdm3 = make_shared<RDM<3>>(norb);
  auto rdm4 = nullptr; // make_shared<RDM<4>>(norb);

  // first make <I|E_ij|bra>
  auto dbra = make_shared<Dvec>(cbra->det(), norb*norb);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // second make <J|E_kl|I><I|E_ij|bra> - delta_li <J|E_kj|bra> = <J|E_kl,ij|bra> = <J|E_ij,kl|bra>
  auto ebra = make_shared<Dvec>(cbra->det(), norb*norb*norb*norb);
  auto tmp = make_shared<Dvec>(cbra->det(), norb*norb);
  int ijkl = 0;
  int ij = 0;
  for (auto iter = dbra->dvec().begin(); iter != dbra->dvec().end(); ++iter, ++ij) {
    const int j = ij/norb;
    const int i = ij-j*norb;
    tmp->zero();
    sigma_2a1(*iter, tmp);
    sigma_2a2(*iter, tmp);
    int kl = 0;
    for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
      *ebra->data(ijkl) = **t;
      const int l = kl/norb;
      const int k = kl-l*norb;
      if (l == i) *ebra->data(ijkl) -= *dbra->data(k+j*norb);
    }
  }

  //ADDED: make dket & eket
  //dket = <I|E_ij|ket>
  shared_ptr<Dvec> dket;
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb*norb);
    dket->zero();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }
  //eket = <I|E_ij,kl|ket>
//auto eket = make_shared<Dvec>(cket->det(), norb*norb*norb*norb);
//auto tmpk = make_shared<Dvec>(cket->det(), norb*norb);
//ijkl = 0;
//ij = 0;
//for (auto iter = dket->dvec().begin(); iter != dket->dvec().end(); ++iter, ++ij) {
//  const int j = ij/norb;
//  const int i = ij-j*norb;
//  tmpk->zero();
//  sigma_2a1(*iter, tmpk);
//  sigma_2a2(*iter, tmpk);
//  int kl = 0;
//  for (auto t = tmpk->dvec().begin(); t != tmpk->dvec().end(); ++t, ++ijkl, ++kl) {
//    *eket->data(ijkl) = **t;
//    const int l = kl/norb;
//    const int k = kl-l*norb;
//    if (l == i) *eket->data(ijkl) -= *dket->data(k+j*norb);
//  }
//}

  //RDM1(not used) & RDM2 build
  auto rdm1 = make_shared<RDM<1>>(norb); //unused
  auto rdm2 = make_shared<RDM<2>>(norb); // <bra|E_ij,kl|ket>
  tie(rdm1,rdm2) = compute_rdm12_last_step(dbra, dket, cbra);
  //END

  // size of the RI space
  const size_t nri = ebra->lena() * ebra->lenb();
  assert(nri == dbra->lena()*dbra->lenb());
  //ADDED
//assert(nri == eket->lena() * eket->lenb() && nri == dket->lena() * dket->lenb());
  //END

  // first form <bra|E_ij,kl|I><I|E_mn|ket>
  {
    auto tmp3 = make_shared<RDM<3>>(norb);
    dgemm_("T", "N", ebra->ij(), dket->ij(), nri, 1.0, ebra->data(), nri, dket->data(), nri, 0.0, tmp3->data(), ebra->ij());
  
    //reorder oprators, note p25
    SMITH::sort_indices<1,0,3,2,4,5, 0,1,1,1>(tmp3->data(), rdm3->data(), norb, norb, norb, norb, norb, norb);

    // then perform Eq. 49 of JCP 89 5803 (Werner's MRCI paper)
    for (int n = 0; n != norb; ++n) 
    for (int l = 0; l != norb; ++l) 
    for (int k = 0; k != norb; ++k) 
    for (int j = 0; j != norb; ++j) 
    for (int i = 0; i != norb; ++i) {
      rdm3->element(i, j, k, l, j, n) -= rdm2->element(i,n,k,l);
      rdm3->element(i, j, k, l, l, n) -= rdm2->element(i,j,k,n);
    }
  }
  cout << "rdm3" << endl; cout.flush();

#if 0
  // 4RDM <0|E_ij,kl|I><I|E_mn,op|0>
  {
    {
      auto tmp4 = make_shared<RDM<4>>(norb);
      dgemm_("T", "N", ebra->ij(), eket->ij(), nri, 1.0, ebra->data(), nri, eket->data(), nri, 0.0, tmp4->data(), ebra->ij());
      SMITH::sort_indices<1,0,3,2,4,5,6,7,0,1,1,1>(tmp4->data(), rdm4->data(), norb, norb, norb, norb, norb, norb, norb, norb);
      for (int l = 0; l != norb; ++l)
        for (int d = 0; d != norb; ++d)
          for (int k = 0; k != norb; ++k)
            for (int c = 0; c != norb; ++c)
              for (int j = 0; j != norb; ++j)
                for (int b = 0; b != norb; ++b)
                  for (int i = 0; i != norb; ++i)
                    for (int a = 0; a != norb; ++a) {
                      if (c == i && d == j) rdm4->element(a,i,b,j,c,k,d,l) -= rdm2->element(a,k,b,l); //TODO need check
                      if (c == j && d == i) rdm4->element(a,i,b,j,c,k,d,l) -= rdm2->element(a,l,b,k); // same
                      if (c == i)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,k,b,j,d,l);
                      if (c == j)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,i,b,k,d,l);
                      if (d == i)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,l,b,j,c,k);
                      if (d == j)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,i,b,l,c,k);
                    }
    }
  }
  cout << "rdm4" << endl; cout.flush();
#endif

/*
  // Checking 4RDM by comparing with 3RDM
  const int nelea = cbra->det()->nelea();
  const int neleb = cbra->det()->neleb();
  auto debug = make_shared<RDM<3>>(*rdm3);
  cout << "printing out rdm" << endl;
  for (int l = 0; l != norb; ++l)
    for (int d = 0; d != norb; ++d)
      for (int k = 0; k != norb; ++k)
        for (int c = 0; c != norb; ++c)
          for (int j = 0; j != norb; ++j)
            for (int b = 0; b != norb; ++b)
    for (int i = 0; i != norb; ++i) {
      debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * rdm4->element(i,i,b,j,c,k,d,l);
//    debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(b,j,i,i,c,k,d,l);
//    debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(b,j,c,k,i,i,d,l);
//    debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(b,j,c,k,d,l,i,i);
    }
  debug->print(1.0e-8);
  cout << "printing out rdm - end" << endl;
*/

  return tie(rdm3, rdm4);
}


//***************************************************************************************************************
std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
ASD_CAS::compute_rdm12_monomer (std::pair<int,int> offset, std::array<Dvec,4>& fourvecs) const {
//returns monomer 1&2RDMs
//cf. ASD_RAS, ASD_DistCAS, ASD_DistRAS for the same function with different vector types
//***************************************************************************************************************
  std::cout << "ASD_CAS: compute_rdm12_monomer called" << std::endl;

  //Dvec decomposition
  auto& A  = fourvecs[0]; // <I'>
  auto& Ap = fourvecs[1]; // |I>

  auto& B  = fourvecs[2]; // <J'|
  auto& Bp = fourvecs[3]; // |J>

  //Offsets
  const int ioff = std::get<0>(offset);
  const int joff = std::get<1>(offset);

  const int nactA = dimer_->active_refs().first->nact(); //dimer_ (ASD_base)
  const int nactB = dimer_->active_refs().second->nact();

  auto rdm1A = std::make_shared<RDM<1>>(nactA);
  auto rdm2A = std::make_shared<RDM<2>>(nactA);
  auto rdm1B = std::make_shared<RDM<1>>(nactB);
  auto rdm2B = std::make_shared<RDM<2>>(nactB);

  rdm1A->zero(); rdm2A->zero();
  rdm1B->zero(); rdm2B->zero();
  const int nstA  = A.ij();
  const int nstAp = Ap.ij();
  std::cout << "<I'| x |I> = " << nstA << " x " << nstAp << std::endl;

  const int nstB  = B.ij();
  const int nstBp = Bp.ij();
  std::cout << "<J'| x |J> = " << nstA << " x " << nstAp << std::endl;

  assert(nstA == nstB && nstAp == nstBp);

  //currently, Diagonal subspace only
  assert(nstA == nstAp);

  //ground state only
  //const int istate = 0;

  
  //Monomer A
//for (int j = 0; j != nstB; ++j) { // <J'|
//  for (int i = 0; i != nstA; ++i) { // <I'|
//    const int ij = i + (j*nstA); //cf. dimerindex()

//    for (int jp = 0; jp != nstBp; ++jp) { // |J>
//      for (int ip = 0; ip != nstAp; ++ip) { // |I>
//        const int ijp = ip + (jp*nstAp);
//        const double coef = adiabats_->element(ioff+ij,0) * adiabats_->element(joff+ijp,0); // C_(I'J') * C_(IJ) TODO: 0 = ground state only

//        if(j == jp) { //delta_J'J
//          shared_ptr<RDM<1>> r1;
//          shared_ptr<RDM<2>> r2;
//          tie(r1,r2) = compute_rdm12_from_civec(A.data(i), Ap.data(ip)); // <I'|E(op)|I>
//          r1->scale(coef);
//          *rdm1A += *r1;
//          r2->scale(coef);
//          *rdm2A += *r2;
//        }

//        if(i == ip) { //delta_I'I
//          shared_ptr<RDM<1>> r1;
//          shared_ptr<RDM<2>> r2;
//          tie(r1,r2) = compute_rdm12_from_civec(B.data(j), Bp.data(jp)); // <J'|E(op)|J>
//          r1->scale(coef);
//          *rdm1B += *r1;
//          r2->scale(coef);
//          *rdm2B += *r2;
//        }

//      } //ip
//    } //jp

//  } //i
//} //j
  
  //NEW
  //MonomerA
  for (int i = 0; i != nstA; ++i) {//<I'|
    for (int ip = 0; ip != nstAp; ++ip) {// |I>
  
      shared_ptr<RDM<1>> r1;
      shared_ptr<RDM<2>> r2;
      tie(r1,r2) = compute_rdm12_from_civec(A.data(i), Ap.data(ip)); // <I'|E(op)|I>

      double csum = 0.0; //coeff sum
      for (int j = 0; j != nstB; ++j) { // delta_J'J
        const int ij  = i  + (j*nstA);
        const int ijp = ip + (j*nstAp);
        csum += adiabats_->element(ioff+ij,0) * adiabats_->element(joff+ijp,0);
      } //I'I

      r1->scale(csum);
      *rdm1A += *r1;
      r2->scale(csum);
      *rdm2A += *r2;
    } //|I>
  } //<I'|
  //MonomerB
  for (int j = 0; j != nstB; ++j) {//<J'|
    for (int jp = 0; jp != nstBp; ++jp) {// |J>
  
      shared_ptr<RDM<1>> r1;
      shared_ptr<RDM<2>> r2;
      tie(r1,r2) = compute_rdm12_from_civec(B.data(j), Bp.data(jp)); // <J'|E(op)|J>

      double csum = 0.0; //coeff sum
      for (int i = 0; i != nstA; ++i) { // delta_I'I
        const int ij  = i + (j*nstA);
        const int ijp = i + (jp*nstAp);
        csum += adiabats_->element(ioff+ij,0) * adiabats_->element(joff+ijp,0);
      } //I'I

      r1->scale(csum);
      *rdm1B += *r1;
      r2->scale(csum);
      *rdm2B += *r2;
    } //|J>
  } //<J'|
  //END NEW

  auto out1 = std::make_shared<RDM<1>>(nactA+nactB);
  out1->zero();
  {
    //Monomer A
    auto low = {0,0};
    auto up  = {nactA,nactA};
    auto outv = make_rwview(out1->range().slice(low,up), out1->storage());
    copy(rdm1A->begin(), rdm1A->end(), outv.begin());
  }
  {
    //Monomer B
    auto low = {nactA,nactA};
    auto up  = {nactA+nactB,nactA+nactB};
    auto outv = make_rwview(out1->range().slice(low,up), out1->storage());
    copy(rdm1B->begin(), rdm1B->end(), outv.begin());
  }
  auto out2 = std::make_shared<RDM<2>>(nactA+nactB);
  out2->zero();
  {
    //Monomer A
    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto outv = make_rwview(out2->range().slice(low,up), out2->storage());
    copy(rdm2A->begin(), rdm2A->end(), outv.begin());
  }
  {
    //Monomer B
    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactA+nactB,nactA+nactB,nactA+nactB,nactA+nactB};
    auto outv = make_rwview(out2->range().slice(low,up), out2->storage());
    copy(rdm2B->begin(), rdm2B->end(), outv.begin());
  }

  return std::make_tuple(out1, out2);

}

//TODO
//below functions are taken from fci/fci_rdm.cc & fci/knowles_compute.cc
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> 
ASD_CAS::compute_rdm12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {

  //ADDED
  const int norb = cbra->det()->norb();
  assert(*cbra->det() == *cket->det());
  //END

  // since we consider here number conserving operators...
  auto dbra = make_shared<Dvec>(cbra->det(), norb*norb);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb*norb);
    dket->zero();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  return compute_rdm12_last_step(dbra, dket, cbra);
}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> 
ASD_CAS::compute_rdm12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {

  //ADDED
  const int norb = cibra->det()->norb();
  //END

  const int nri = dbra->lena()*dbra->lenb();
  const int ij  = norb*norb;

  if (nri != dket->lena()*dket->lenb())
    throw logic_error("FCI::compute_rdm12_last_step called with inconsistent RI spaces");

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  auto rdm1 = make_shared<RDM<1>>(norb);
  dgemv_("T", nri, ij, 1.0, dket->data(0)->data(), nri, cibra->data(), 1, 0.0, rdm1->data(), 1);
  // 2RDM
  // \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm2 = make_shared<RDM<2>>(norb);
  dgemm_("T", "N", ij, ij, nri, 1.0, dbra->data(0)->data(), nri, dket->data(0)->data(), nri, 0.0, rdm2->data(), ij);

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb**2).
  unique_ptr<double[]> buf(new double[norb*norb]);
  for (int i = 0; i != norb; ++i) {
    for (int k = 0; k != norb; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb*norb, buf.get());
      blas::transpose(buf.get(), norb, norb, &rdm2->element(0,0,k,i));
    }
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb; ++i)
    for (int k = 0; k != norb; ++k)
      for (int j = 0; j != norb; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}

void ASD_CAS::sigma_2a1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {
  assert(d->det() == cc->det());
  const int lb = d->lenb();
  const int ij = d->ij();
  const double* const source_base = cc->data();
  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = d->data(ip)->data();
    for (auto& iter : cc->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = target_base + iter.source*lb;
      daxpy_(lb, sign, source_base + iter.target*lb, 1, target_array, 1);
    }
  }
}

void ASD_CAS::sigma_2a2(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {
  assert(d->det() == cc->det());
  const int la = d->lena();
  const int ij = d->ij();
  for (int i = 0; i < la; ++i) {
    const double* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      double* const target_array0 = d->data(ip)->element_ptr(0, i);
      for (auto& iter : cc->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}
