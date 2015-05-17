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

#include <src/util/prim_op.h>
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


tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_monomer(shared_ptr<const RASDvec> civec, const int i) const {
  shared_ptr<const RASCivec> cbra = civec->data(i);
  const int norb = cbra->det()->norb();

  Timer mtime;

  auto dbra = make_shared<RASDvec>(cbra->det(), norb*norb);
  dbra->zero();
  sigma_2a(cbra, dbra);

  //cbra == cket
  shared_ptr<RASDvec> dket;
  dket = dbra;

  std::cout << "  o single monomer RDM - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;

  return compute_rdm12_last_step(cbra, dbra);
}

void ASD_RAS::sigma_2a(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& block : cc->blocks()) {
    if (block) {
      const size_t a_offset = block->stringsa()->offset();
      const size_t b_offset = block->stringsb()->offset();

      for (size_t ia = 0, ab = 0; ia < block->stringsa()->size(); ++ia) {
        const auto abit = block->string_bits_a(ia);
        for (size_t jb = 0; jb < block->stringsb()->size(); ++jb, ++ab) {
          const auto bbit = block->string_bits_b(jb);
          double coef = block->element(ab);

          for (auto& phi : det->phib(jb + b_offset)) {
            assert(phi.target == jb + b_offset);
            const auto sbit = det->string_bits_b(phi.source);
            if(det->allowed(abit,sbit))
              d->data(phi.ij)->element(sbit,abit) += static_cast<double>(phi.sign) * coef;
          }

          for (auto& phi : det->phia(ia + a_offset)) {
            assert(phi.target == ia + a_offset);
            const auto sbit = det->string_bits_a(phi.source);
            if(det->allowed(sbit,bbit))
              d->data(phi.ij)->element(bbit,sbit) += static_cast<double>(phi.sign) * coef;
          }

        }
      }
    }
  }

}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_RAS::compute_rdm12_last_step(shared_ptr<const RASCivec> cibra, shared_ptr<const RASDvec> dbra) const {

  Timer mtime;
  const int norb = cibra->det()->norb();

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  shared_ptr<const RASDeterminants> det = cibra->det();
  auto rdm1 = make_shared<RDM<1>>(norb);
  {
    double* target = rdm1->element_ptr(0,0);
    for (int ij = 0; ij != norb*norb; ++ij, ++target){
      *target = cibra->dot_product(dbra->data(ij));
    }
    std::cout << "      o 1RDM - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;
  }

  auto rdm2 = make_shared<RDM<2>>(norb);
  {
    auto rdmt = rdm2->clone();
    double* target = rdmt->element_ptr(0,0,0,0);
    for (int kl = 0; kl != norb*norb; ++kl) // E_kl |ket>
      for (int ij = 0; ij != norb*norb; ++ij, ++target) // <bra| E_ji = [E_ij |bra>]^+
        *target = dbra->data(ij)->dot_product(dbra->data(kl)); // (ij,kl) has E_ji*E_kl

    std::cout << "      o 2RDM - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;

    unique_ptr<double[]> buf(new double[norb*norb]);
    for (int i = 0; i != norb; ++i)
      for (int k = 0; k != norb; ++k) {
        copy_n(&rdmt->element(0,0,k,i), norb*norb, buf.get());
        blas::transpose(buf.get(), norb, norb, &rdmt->element(0,0,k,i));
      }
    std::cout << "      o tran - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;

    //Swap kl,ij -> ij,kl
    sort_indices<2,3,0,1, 0,1, 1,1>(rdmt->data(), rdm2->data(), norb, norb, norb, norb);
    std::cout << "      o sort - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;

    // put in diagonal into 2RDM
    // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
    for (int i = 0; i != norb; ++i)
      for (int k = 0; k != norb; ++k)
        for (int j = 0; j != norb; ++j)
          rdm2->element(j,k,k,i) -= rdm1->element(j,i);
    std::cout << "      o diag - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;

    //RDM2 symmetrize (out-of-excitation free parts are copied)
    for (int i = 0, ij = 0; i != norb; ++i)
      for (int j = 0; j != norb; ++j, ++ij)
          for (int k = 0, kl = 0; k != norb; ++k)
            for (int l = 0; l != norb; ++l, ++kl)
              if (kl > ij) rdm2->element(i,j,k,l) = rdm2->element(k,l,i,j);
    std::cout << "      o sym - " << std::setw(9) << std::fixed << std::setprecision(5) << mtime.tick() << std::endl;
  }

  return tie(rdm1, rdm2);
}

shared_ptr<RASDvec> ASD_RAS::contract_I(shared_ptr<const RASDvec> A, shared_ptr<Matrix> adiabats, int ioff, int nstA, int nstB, int kst) const {
  auto out = make_shared<RASDvec>(A->det(), nstB);
  for (int ij = 0; ij != nstB; ++ij) {
    out->data(ij)->zero();
  }

  for (int j = 0; j != nstB; ++j) {
    for (int i = 0; i != nstA; ++i) {
      const int ij  = i  + (j*nstA);
      double u_ij = adiabats->element(ioff+ij,kst);

      out->data(j)->ax_plus_y(u_ij, A->data(i));

    }
  }
  return out;
}

shared_ptr<RASDvec> ASD_RAS::contract_J(shared_ptr<const RASDvec> B, shared_ptr<Matrix> adiabats, int ioff, int nstA, int nstB, int kst) const {
  auto out = make_shared<RASDvec>(B->det(), nstA);
  for (int ij = 0; ij != nstA; ++ij) {
    out->data(ij)->zero();
  }

  for (int i = 0; i != nstA; ++i) {
    for (int j = 0; j != nstB; ++j) {
      const int ij  = i  + (j*nstA);
      double u_ij = adiabats->element(ioff+ij,kst);

      out->data(i)->ax_plus_y(u_ij, B->data(j));

    }
  }
  return out;
}
