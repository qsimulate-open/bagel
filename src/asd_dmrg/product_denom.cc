//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/product_denom.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/asd_dmrg/product_rasci.h>
#include <src/math/btas_interface.h>
#include <src/ras/denomtask.h>

using namespace std;
using namespace bagel;

void ProductRASCI::construct_denom() {
  Timer denom_t;

  // allocate denom_
  denom_ = make_shared<ProductRASCivec>(space_, left_, nelea_, neleb_);

  const int lnorb = left_->norb();
  const int rnorb = norb_;

  // first compute pure Block terms
  map<BlockKey, shared_ptr<VectorB>> pure_block;
  for (auto& b : left_->blocks()) {
    auto out = make_shared<VectorB>(b.nstates);
    shared_ptr<const Matrix> h2e = left_->h2e(b.key());
    const int nstates = b.nstates;
    for (int i = 0; i < nstates; ++i)
      (*out)(i) = h2e->element(i,i);
    pure_block.emplace(b.key(), out);

    const Matrix block1e(*jop_->monomer_jop<1>()->mo1e()->matrix());

    // alpha-alpha part
    const btas::TensorView4<double> alphaview(btas::make_view(btas::CRange<4>(nstates, nstates, lnorb, lnorb), left_->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha}).at({b.key(), b.key()}).data->storage()));
    const btas::TensorView4<double> betaview(btas::make_view(btas::CRange<4>(nstates, nstates, lnorb, lnorb), left_->coupling({GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta}).at({b.key(), b.key()}).data->storage()));
    for (int ist = 0; ist < nstates; ++ist) {
      for (int i = 0; i < lnorb; ++i) {
        for (int j = 0; j < lnorb; ++j) {
          (*out)(ist) += (alphaview(ist,ist,i,j)+betaview(ist,ist,i,j)) * block1e(i,j);
        }
      }
    }
  }
  denom_t.tick_print("denom: pure block");

  // now compute pure RAS terms
  map<BlockKey, shared_ptr<RASCivec>> pure_ras;
  {
    VectorB h(rnorb);
    Matrix jop(rnorb, rnorb);
    Matrix kop(rnorb, rnorb);

    shared_ptr<const MOFile> rasjop = jop_->monomer_jop<0>();
    for (int i = 0; i != rnorb; ++i) {
      for (int j = 0; j <= i; ++j) {
        jop(i,j) = jop(j,i) = 0.5*rasjop->mo2e_hz(j, i, j, i);
        kop(i,j) = kop(j,i) = 0.5*rasjop->mo2e_hz(j, i, i, j);
      }
      h(i) = rasjop->mo1e(i,i);
    }

    TaskQueue<RAS::DenomTask> tasks;
    for (auto& sec : denom_->sectors()) {
      shared_ptr<const RASDeterminants> det = sec.second->det();
      auto civec = make_shared<RASCivec>(det);
      for (auto& iblock : civec->blocks()) {
        if (!iblock) continue;
        double* iter = iblock->data();
        for (auto& ia : *iblock->stringsa()) {
          tasks.emplace_back(iter, ia, iblock->stringsb(), jop.data(), kop.data(), h.data());
          iter += iblock->lenb();
        }
      }
      pure_ras.emplace(sec.first, civec);
    }
    tasks.compute();
  }
  denom_t.tick_print("denom: pure ras");

  // finally, pull it all together and compute the mixed terms
  {
    btas::Tensor3<double> jop(rnorb, lnorb, lnorb);
    btas::Tensor3<double> kop(rnorb, lnorb, lnorb);
    {
      const btas::TensorView4<double> fulljop(btas::make_view(btas::CRange<4>(rnorb,rnorb,lnorb,lnorb), jop_->coulomb_matrix<1,0,1,0>()->storage()));
      const btas::TensorView4<double> fullkop(btas::make_view(btas::CRange<4>(rnorb,rnorb,lnorb,lnorb), jop_->coulomb_matrix<1,1,0,0>()->storage()));
      for (int i = 0; i < lnorb; ++i) {
        for (int j = 0; j < lnorb; ++j) {
          for (int k = 0; k < rnorb; ++k) {
            jop(k,i,j) = fulljop(i,j,k,k);
            kop(k,i,j) = fullkop(i,j,k,k);
          }
        }
      }
    }

    // TODO: This should be threaded
    for (auto& sec: denom_->sectors()) {
      const int nstates = sec.second->nstates();
      vector<RASCivecView> vecs = sec.second->civecs();

      const btas::TensorView4<double> alphaview(btas::make_view(btas::CRange<4>(nstates, nstates, lnorb, lnorb), left_->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha}).at({sec.first, sec.first}).data->storage()));
      const btas::TensorView4<double> betaview(btas::make_view(btas::CRange<4>(nstates, nstates, lnorb, lnorb), left_->coupling({GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta}).at({sec.first, sec.first}).data->storage()));

      for (int i = 0; i < nstates; ++i) {
        // initialize to sum of pure terms
        RASCivecView civec = vecs[i];
        civec.fill((*pure_block[sec.first])(i));
        civec.ax_plus_y(1.0, *pure_ras[sec.first]);

        for (auto& block : civec.blocks()) {
          if (!block) continue;
          const size_t la = block->lena();
          const size_t lb = block->lenb();
          for (size_t ia = 0; ia < la; ++ia) {
            bitset<nbit__> abit = block->string_bits_a(ia);
            double* const data_base = block->data() + lb*ia;
            for (int p = 0; p < lnorb; ++p) {
              for (int q = 0; q < lnorb; ++q) {
                double ja_pq = 0.0; // ja_pq = \sum_r n_r^alpha (pq|rr)
                double ka_pq = 0.0; // ka_pq = \sum_r n_r^alpha (pr|qr)
                for (int r = 0; r < rnorb; ++r) {
                  ja_pq += jop(r,p,q)*abit[r];
                  ka_pq += kop(r,p,q)*abit[r];
                }

                for (size_t ib = 0; ib < lb; ++ib) {
                  bitset<nbit__> bbit = block->string_bits_b(ib);

                  double j_pq = ja_pq; // j_pq = \sum_r n_r^beta (pq|rr) + ja_pq
                  double kb_pq = 0.0; // kb_pq = \sum_r n_r^beta (pr|qr)
                  for (int r = 0; r < rnorb; ++r) {
                    j_pq += jop(r,p,q)*bbit[r];
                    kb_pq += kop(r,p,q)*bbit[r];
                  }
                  data_base[ib] += alphaview(i,i,p,q) * (j_pq - ka_pq) + betaview(i,i,p,q) * (j_pq - kb_pq);
                }
              }
            }
          }
        }
      }
    }
  }
  denom_t.tick_print("denom: mixed");
}
