//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_ras_sigma.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/util/prim_op.h>
#include <src/asd/asd_ras.h>
#include <src/ci/ras/form_sigma.h>

using namespace std;
using namespace bagel;

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

  auto dbra = make_shared<RASDvec>(cbra->det(), norb*norb);
  dbra->zero();
  sigma_2a(cbra, dbra);

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
          const double coef = block->element(ab);

          if (fabs(coef) < numerical_zero__) continue;
          for (auto& phi : det->phib(jb + b_offset)) {
            assert(phi.target == jb + b_offset);
            const auto sbit = det->string_bits_b(phi.source);
            if(det->allowed(abit,sbit))
              d->data(phi.ij)->element(sbit,abit) += static_cast<double>(phi.sign) * coef; // i^+ j => (j,i)
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
  const int norb = cibra->det()->norb();
  const int nri = dbra->data(0)->size();
  auto cimat = make_shared<Matrix>(nri,norb*norb);
  for (int ij = 0; ij != norb*norb; ++ij)
    copy_n(dbra->data(ij)->data(), nri, cimat->element_ptr(0,ij));

  // 1RDM
  auto rdm1 = make_shared<RDM<1>>(norb);
  dgemv_("T", nri, norb*norb, 1.0, cimat->element_ptr(0,0), nri, cibra->data(), 1, 0.0, rdm1->data(), 1);

  // 2RDM
  auto rdm2 = make_shared<RDM<2>>(norb);
  {
    auto rdmt = rdm2->clone();
    dgemm_("T", "N", norb*norb, norb*norb, nri, 1.0, cimat->element_ptr(0,0), nri, cimat->element_ptr(0,0), nri, 0.0, rdmt->data(), norb*norb);

    unique_ptr<double[]> buf(new double[norb*norb]);
    for (int i = 0; i != norb; ++i)
      for (int k = 0; k != norb; ++k) {
        copy_n(rdmt->element_ptr(0,0,k,i), norb*norb, buf.get());
        blas::transpose(buf.get(), norb, norb, rdmt->element_ptr(0,0,k,i));
      }

    sort_indices<2,3,0,1, 0,1, 1,1>(rdmt->data(), rdm2->data(), norb, norb, norb, norb);

    // put in diagonal into 2RDM
    // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
    for (int i = 0; i != norb; ++i)
      for (int k = 0; k != norb; ++k)
        for (int j = 0; j != norb; ++j)
          rdm2->element(j,k,k,i) -= rdm1->element(j,i);

    //RDM2 symmetrize (out-of-excitation free parts are copied)
    for (int i = 0, ij = 0; i != norb; ++i)
      for (int j = 0; j != norb; ++j, ++ij)
          for (int k = 0, kl = 0; k != norb; ++k)
            for (int l = 0; l != norb; ++l, ++kl)
              if (kl > ij) rdm2->element(i,j,k,l) = rdm2->element(k,l,i,j);
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
