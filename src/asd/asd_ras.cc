//
// BAGEL - Parallel electron correlation program.
// Filename: asd_ras_sigma.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/asd/asd_ras.h>
#include <src/ci/ras/form_sigma.h>

using namespace std;
using namespace bagel;

ASD_RAS::ASD_RAS(const shared_ptr<const PTree> input, shared_ptr<Dimer> dimer, shared_ptr<DimerRAS> cispace)
  : ASD<RASDvec>(input, dimer, cispace) {

}


shared_ptr<RASDvec> ASD_RAS::form_sigma(shared_ptr<const RASDvec> ccvec, shared_ptr<const MOFile> jop) const {
  FormSigmaRAS form;
  vector<int> conv(ccvec->ij(), static_cast<int>(false));
  return form(ccvec, jop, conv);
}

shared_ptr<RASDvec> ASD_RAS::form_sigma_1e(shared_ptr<const RASDvec> ccvec, const double* modata) const {
  FormSigmaRAS form;
  const int norb = ccvec->det()->norb();
  auto mo1e = make_shared<Matrix>(norb, norb);
  copy_n(modata, norb*norb, mo1e->data());
  return form(ccvec, mo1e, nullptr, vector<int>(ccvec->ij(), static_cast<int>(false)));
}

std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>, std::shared_ptr<RDM<4>>> ASD_RAS::compute_rdm34_monomer (std::pair<int,int> offset, std::array<RASDvec,4>& fourvecs) const {
  std::cout << "ASD_RAS: compute_rdm34_monomer called" << std::endl;
  assert(false);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_monomer (pair<int,int> offset, array<RASDvec,4>& fourvecs) const {
  cout << "ASD_RAS:: compute_rdm12_monomer called" << endl;
  //based on asd_cas.cc(compute_rdm12_monomer)
  //TODO: template this

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
ASD_RAS::compute_rdm12_from_civec(shared_ptr<const RASCivec> cbra, shared_ptr<const RASCivec> cket) const {

  //ADDED
  const int norb = cbra->det()->norb();
  assert(*cbra->det() == *cket->det());
  //END

  // since we consider here number conserving operators...
  auto dbra = make_shared<RASDvec>(cbra->det(), norb*norb);
//dbra->zero();
//sigma_2a1(cbra, dbra);
//sigma_2a2(cbra, dbra);

  shared_ptr<RASDvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<RASDvec>(cket->det(), norb*norb);
//  dket->zero();
//  sigma_2a1(cket, dket);
//  sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

//return compute_rdm12_last_step(dbra, dket, cbra);
  auto rdm1 = make_shared<RDM<1>>(norb);
  auto rdm2 = make_shared<RDM<2>>(norb);
  return tie(rdm1,rdm2);
}

#if 0

void ASD_RAS::sigma_2a1(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
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

void ASD_RAS::sigma_2a2(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
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

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> 
ASD_CAS::compute_rdm12_last_step(shared_ptr<const RASDvec> dbra, shared_ptr<const RASDvec> dket, shared_ptr<const RASCivec> cibra) const {

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

#endif
