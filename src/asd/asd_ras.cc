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


tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_monomer(shared_ptr<const RASDvec> civec, const int i, shared_ptr<const RASDvec> cipvec,  const int ip) const {
  cout << "compute_rdm12_from_civec" << endl;
  shared_ptr<const RASCivec> cbra = civec->data(i);
  shared_ptr<const RASCivec> cket = cipvec->data(ip);
  const int norb = cbra->det()->norb();
  assert(*cbra->det() == *cket->det());

  auto dket = make_shared<RASDvec>(cket->det(), norb*norb);
  for (int ij = 0; ij != norb*norb; ++ij) {
    dket->data(ij)->zero();
  }

  sigma_2a1(cket, dket);
  sigma_2a2(cket, dket);

  auto eket = make_shared<RASDvec>(cket->det(), norb*norb*norb*norb);
  for (int ij = 0; ij != norb*norb*norb*norb; ++ij) {
    eket->data(ij)->zero();
  }

  sigma_2a1_aa(cket, eket); //aa
  sigma_2a2_bb(cket, eket); //bb
  sigma_2a3_ba(cket, eket); //ba TODO *2 to remove ab?
  sigma_2a4_ab(cket, eket); //ab

  return compute_rdm12_last_step(cbra, dket, eket);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_last_step(shared_ptr<const RASCivec> cibra, shared_ptr<const RASDvec> dket, shared_ptr<const RASDvec> eket) const {

  const int norb = cibra->det()->norb();
  cout << "last-step entered.." << endl;

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  shared_ptr<const RASDeterminants> det = cibra->det();
  auto rdm1 = make_shared<RDM<1>>(norb);
  for (int i = 0, ij = 0; i != norb; ++i){
    for (int j = 0; j!= norb; ++j, ++ij) {

      for (auto& cblock : cibra->blocks()) {
        if (!cblock) continue;
  
        for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
          for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
            auto abit = cblock->stringsa()->strings(ca);
            auto bbit = cblock->stringsb()->strings(cb);
            rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cibra->element(bbit,abit); 
          }
        }
      }
    }
  }
  cout << "1rdm done.." << endl;

  auto rdm2 = make_shared<RDM<2>>(norb);
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
                assert(det->allowed(abit,bbit));
                rdm2->element(k,l,i,j) += cibra->element(bbit,abit) * eket->data(ijkl)->element(bbit,abit);
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
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  cout << "2rdm done.." << endl;

  //RDM2 symmetrize (out-of-excitation free parts are copied)
  for (int i = 0, ij = 0; i != norb; ++i) 
    for (int j = 0; j != norb; ++j, ++ij) 
      for (int k = 0, kl = 0; k != norb; ++k)
        for (int l = 0; l != norb; ++l, ++kl)
          if (kl > ij) rdm2->element(i,j,k,l) = rdm2->element(k,l,i,j);

  return tie(rdm1, rdm2);
}



void ASD_RAS::sigma_2a1(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspaceb()) { 
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const auto bbit = ispace->strings(ib);

      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) {

          for (auto& phi : det->phia(ja+offset)) {
            assert(phi.target == ja+offset);
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
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspacea()) { // alpha determinant space
    for (size_t ia = 0; ia != ispace->size(); ++ia) { // determinants associated with a given space
      const auto abit = ispace->strings(ia);

      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) { // determinants associated with a given space

          for (auto& phi : det->phib(jb+offset)) {
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_b(phi.source);
            const auto tbit = det->string_bits_b(phi.target);
            const auto ij = phi.ij;
            if(!det->allowed(abit,tbit)) continue;
            if(!det->allowed(abit,sbit)) continue;
            d->data(ij)->element(sbit,abit) += sign * cc->element(tbit,abit);
          }
        }
      }
    }
  }

}


void ASD_RAS::sigma_2a1_aa(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
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
            const auto tbit = det->string_bits_a(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(tbit,bbit)) continue;

            for (auto& phi2 : det->phia(phi.source)) {
              assert(phi2.target == phi.source);
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


void ASD_RAS::sigma_2a2_bb(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = cc->det()->norb();

  for (auto& ispace : *det->stringspacea()) { 
    for (size_t ia = 0; ia != ispace->size(); ++ia) {
      const auto abit = ispace->strings(ia);

      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset(); //scan through beta-det
        for (size_t jb = 0; jb != jspace->size(); ++jb) { 

          for (auto& phi : det->phib(jb+offset)) { //first E
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
          //const auto ibit = det->string_bits_b(phi.source); //intermediate
            const auto tbit = det->string_bits_b(phi.target); //coeff
            const auto kl = phi.ij;
            if(!det->allowed(abit,tbit)) continue; //coeff

            for (auto& phi2 : det->phib(phi.source)) {
              assert(phi2.target == phi.source);
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


void ASD_RAS::sigma_2a3_ba(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = cc->det()->norb();

  for (auto& ispace : *det->stringspaceb()) { 
    const size_t boffset = ispace->offset();
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const auto bbit = ispace->strings(ib); //bitset<nbit__>

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
              d->data(ij + kl*norb*norb)->element(sbbit,sbit) += sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

}


void ASD_RAS::sigma_2a4_ab(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = cc->det()->norb();

  for (auto& ispace : *det->stringspacea()) { 
    const size_t aoffset = ispace->offset();
    for (size_t ia = 0; ia != ispace->size(); ++ia) {
      const bitset<nbit__> abit = ispace->strings(ia);

      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) {  //fix b-string

          for (auto& phi : det->phib(jb+offset)) {
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_b(phi.source); //sbit(beta)
            const auto tbit = det->string_bits_b(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(abit,tbit)) continue; //coeff

            for (auto& phi2 : det->phia(ia+aoffset)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sabit = det->string_bits_a(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(sabit,sbit)) continue; //double replacement ket
              d->data(ij + kl*norb*norb)->element(sbit,sabit) += sign * sign2 * cc->element(tbit,abit);
            }
          }
        }
      }
    }
  }

}
