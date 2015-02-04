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

tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>, shared_ptr<RDM<4>>> ASD_RAS::compute_rdm34_monomer(pair<int,int> offset, array<RASDvec,4>& fourvecs) const {
  cout << "ASD_RAS: compute_rdm34_monomer called" << endl;
  assert(false);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_monomer(pair<int,int> offset, array<RASDvec,4>& fourvecs) const {
  cout << "ASD_RAS:: compute_rdm12_monomer called" << endl;
  //based on asd_cas.cc(compute_rdm12_monomer)
  //TODO: template this

  //Dvec decomposition
  auto& A  = fourvecs[0]; // <I'>
  auto& Ap = fourvecs[1]; // |I>

  auto& B  = fourvecs[2]; // <J'|
  auto& Bp = fourvecs[3]; // |J>

  //Offsets
  const int ioff = get<0>(offset);
  const int joff = get<1>(offset);

  const int nactA = dimer_->active_refs().first->nact(); //dimer_ (ASD_base)
  const int nactB = dimer_->active_refs().second->nact();

  cout << "Active : " << nactA << " " << nactB << endl;

  auto rdm1A = make_shared<RDM<1>>(nactA);
  auto rdm2A = make_shared<RDM<2>>(nactA);
  auto rdm1B = make_shared<RDM<1>>(nactB);
  auto rdm2B = make_shared<RDM<2>>(nactB);

  rdm1A->zero(); rdm2A->zero();
  rdm1B->zero(); rdm2B->zero();
  const int nstA  = A.ij();
  const int nstAp = Ap.ij();
  cout << "<I'| x |I> = " << nstA << " x " << nstAp << endl;

  const int nstB  = B.ij();
  const int nstBp = Bp.ij();
  cout << "<J'| x |J> = " << nstA << " x " << nstAp << endl;

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

  auto out1 = make_shared<RDM<1>>(nactA+nactB);
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
  auto out2 = make_shared<RDM<2>>(nactA+nactB);
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

  return make_tuple(out1, out2);

}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_from_civec(shared_ptr<const RASCivec> cbra, shared_ptr<const RASCivec> cket) const {

  cout << "compute_rdm12_from_civec" << endl;
  const int norb = cbra->det()->norb();
  assert(*cbra->det() == *cket->det());

  // since we consider here number conserving operators...
  auto dbra = make_shared<RASDvec>(cbra->det(), norb*norb);
  for (int ij = 0; ij != norb*norb; ++ij) {
    dbra->data(ij)->zero();
  }

//dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<RASDvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
//if (cbra != cket) {
    dket = make_shared<RASDvec>(cket->det(), norb*norb);
    for (int ij = 0; ij != norb*norb; ++ij) {
      dket->data(ij)->zero();
    }
  //dket->zero();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
//} else {
//  dket = dbra;
//}

  //new
  auto eket = make_shared<RASDvec>(cket->det(), norb*norb*norb*norb);
  for (int ij = 0; ij != norb*norb*norb*norb; ++ij)
    eket->data(ij)->zero();
  sigma_2a1_new(cket, eket); //aa
  sigma_2a2_new(cket, eket); //bb
  sigma_2a3_new(cket, eket); //2*ba

  return compute_rdm12_last_step(dbra, dket, cbra, eket);
}


void ASD_RAS::sigma_2a1(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a1" << endl;
  assert(d->det() == cc->det());
  
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspaceb()) { 
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const bitset<nbit__> bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) { 
          for (auto& phi : det->phia(ja+offset)) {
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_a(phi.source);
            const auto tbit = det->string_bits_a(phi.target);
            const auto ij = phi.ij;
            if(!det->allowed(tbit,bbit)) continue;
            if(!det->allowed(sbit,bbit)) continue;
            d->data(ij)->element(bbit,sbit) += sign * cc->element(bbit,tbit);
          }
        }
      }
    }
  }

}

void ASD_RAS::sigma_2a2(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a2" << endl;
  assert(d->det() == cc->det());
  
//const int nij = norb_*norb_;
  shared_ptr<const RASDeterminants> det = cc->det();
//const size_t lb = det->lenb();

  for (auto& ispace : *det->stringspacea()) { // alpha determinant space
  //const size_t offset = ispace->offset();
    for (size_t ia = 0; ia != ispace->size(); ++ia) { // determinants associated with a given space
      const bitset<nbit__> abit = ispace->strings(ia);

    //for (auto& phi : det->uncompressed_phib(ia+offset)) {
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) { // determinants associated with a given space

          for (auto& phi : det->phib(jb+offset)) {
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_b(phi.source);
            const auto tbit = det->string_bits_b(phi.target);
            const auto ij = phi.ij;
          //if(i == j) {
          //  cout << "diag[" << i << "] : " << phi.source << " " << phi.target << endl;
          //}
            if(!det->allowed(abit,tbit)) continue;
            if(!det->allowed(abit,sbit)) continue;
            d->data(ij)->element(sbit,abit) += sign * cc->element(tbit,abit);
          }
        }
      }
    }
  }
}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  ASD_RAS::compute_rdm12_last_step(shared_ptr<const RASDvec> dbra, shared_ptr<const RASDvec> dket, shared_ptr<const RASCivec> cibra, shared_ptr<const RASDvec> eket) const {

  const int norb = cibra->det()->norb();
  cout << "last-step entered.." << endl;
/*
  const int nri = dbra->lena()*dbra->lenb();
  const int ij  = norb_*norb_;

  if (nri != dket->lena()*dket->lenb())
    throw logic_error("FCI::compute_rdm12_last_step called with inconsistent RI spaces");
*/

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  shared_ptr<const RASDeterminants> det = cibra->det();
  auto rdm1 = make_shared<RDM<1>>(norb);
  rdm1->zero();
//dgemv_("T", nri, ij, 1.0, dket->data(0)->data(), nri, cibra->data(), 1, 0.0, rdm1->data(), 1);
  for (int i = 0, ij = 0; i != norb; ++i){
    for (int j = 0; j!= norb; ++j, ++ij) {

      for (auto& cblock : cibra->blocks()) {
        if (!cblock) continue;
  
        for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
          for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
          //auto cptr = cblock->data() + cab;
          //cout << "[" << ca << "," << cb << "] = " << *cptr << " " <<cblock->element(cab) << endl;
            auto abit = cblock->stringsa()->strings(ca);
            auto bbit = cblock->stringsb()->strings(cb);
          //rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cblock->element(cab);
          //rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cblock->element(bbit,abit); 
            rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cibra->element(bbit,abit); 
          //cout << endl;
          //cout << "check: " << cblock->element(cab) << endl;
          //cout << "       " << cblock->element(bbit,abit) << endl;
          //cout << "       " << cibra->element(bbit,abit) << endl;
          //cout << endl;
          }
        }
      }
    }
  }
  cout << "1rdm done.." << endl;

  auto new2 = make_shared<RDM<2>>(norb);
  {
    //NEW 2RDM
    new2->zero();
    for (int i = 0, ij = 0, ijkl = 0; i != norb; ++i){
      for (int j = 0; j != norb; ++j, ++ij) {
        //fixed ij
        for (int k = 0, kl = 0; k != norb; ++k){
          for (int l = 0; l != norb; ++l, ++kl, ++ijkl) {
            //fixed kl
            
            for (auto& cblock : cibra->blocks()) { //just run over allowed blocks/ TODO: replace this
              if (!cblock) continue;
      
              for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
                for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
                  auto abit = cblock->stringsa()->strings(ca);
                  auto bbit = cblock->stringsb()->strings(cb);
                  new2->element(i,j,k,l) += cibra->element(bbit,abit) * eket->data(ijkl)->element(bbit,abit);
                }
              }
 
            }
          }
        }
 
      }
    }
 
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb; ++i)
    for (int k = 0; k != norb; ++k)
      for (int j = 0; j != norb; ++j)
        new2->element(j,k,k,i) -= rdm1->element(j,i);

  cout << "2rdm done.." << endl;

//return tie(rdm1, rdm2);
  return tie(rdm1, new2);
}

void ASD_RAS::sigma_2a1_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a1_new" << endl;
  //based on sigma_2a2
  
  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = cc->det()->norb();

  for (auto& ispace : *det->stringspaceb()) { 
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const bitset<nbit__> bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) { 
          for (auto& phi : det->phia(ja+offset)) {
            assert(phi.target == ja+offset);
            const double sign = static_cast<double>(phi.sign);
          //const auto ibit = det->string_bits_b(phi.source); //intermediate
            const auto tbit = det->string_bits_a(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(tbit,bbit)) continue;
            for (auto& phi2 : det->phia(phi.source)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbit = det->string_bits_a(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(sbit,bbit)) continue;
              d->data(ij + kl*norb*norb)->element(bbit,sbit) += sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

}

void ASD_RAS::sigma_2a2_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a2_new" << endl;
  
  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = cc->det()->norb();

  for (auto& ispace : *det->stringspacea()) { 
    for (size_t ia = 0; ia != ispace->size(); ++ia) {
      const bitset<nbit__> abit = ispace->strings(ia);
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) { 
          for (auto& phi : det->phib(jb+offset)) {
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
          //const auto ibit = det->string_bits_b(phi.source); //intermediate
            const auto tbit = det->string_bits_b(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(abit,tbit)) continue;
            for (auto& phi2 : det->phib(phi.source)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbit = det->string_bits_b(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(abit,sbit)) continue;
              d->data(ij + kl*norb*norb)->element(sbit,abit) += sign * sign2 * cc->element(tbit,abit);
            }
          }
        }
      }
    }
  }

}

void ASD_RAS::sigma_2a3_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a3_new" << endl;
  //based on sigma_2a2
  
  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = cc->det()->norb();

  for (auto& ispace : *det->stringspaceb()) { 
    const size_t boffset = ispace->offset();
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const bitset<nbit__> bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) {  //fix a-string
          for (auto& phi : det->phia(ja+offset)) {
            assert(phi.target == ja+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_a(phi.source); //sbit(alpha)
            const auto tbit = det->string_bits_a(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(tbit,bbit)) continue; //coeff
            for (auto& phi2 : det->phib(ib+boffset)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbbit = det->string_bits_b(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(sbit,sbbit)) continue; //double replacement ket
            //cout << cc->element(bbit,tbit) << endl;
            //cout << d->data(ij + kl*norb_*norb_)->element(sbbit,sbit) << endl;
            //cout << endl;
              d->data(ij + kl*norb*norb)->element(sbbit,sbit) += 2.0 * sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

}

