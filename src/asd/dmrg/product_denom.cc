//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/product_denom.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/asd/dmrg/product_rasci.h>
#include <src/util/math/btas_interface.h>
#include <src/ci/ras/denomtask.h>

using namespace std;
using namespace bagel;

void ProductRASCI::construct_denom() {
  Timer denom_t;

  // allocate denom_
  denom_ = make_shared<ProductRASCivec>(space_, left_, nelea_, neleb_);

  const int rnorb = norb_;

  // first compute pure Block terms
  map<BlockKey, shared_ptr<VectorB>> pure_block;
  for (auto& b : left_->blocks()) {
    auto out = make_shared<VectorB>(b.nstates);
    shared_ptr<const Matrix> ham = blockops_->ham(b.key());
    const int nstates = b.nstates;
    for (int i = 0; i < nstates; ++i)
      (*out)(i) = ham->element(i,i);
    pure_block.emplace(b.key(), out);
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
    for (auto& sec: denom_->sectors()) {
      const int nstates = sec.second->nstates();
      vector<RASCivecView> vecs = sec.second->civecs();

      Matrix Qaa(rnorb, nstates);
      Matrix Qbb(rnorb, nstates);

      for (int r = 0; r < rnorb; ++r) {
        shared_ptr<const VectorB> qaa_r = blockops_->Q_aa(sec.first,r,r)->diagonal();
        shared_ptr<const VectorB> qbb_r = blockops_->Q_bb(sec.first,r,r)->diagonal();
        for (int i = 0; i < nstates; ++i) {
          Qaa(r, i) = (*qaa_r)(i);
          Qbb(r, i) = (*qbb_r)(i);
        }
      }

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

            double alphaE = 0.0;
            for (int r = 0; r < rnorb; ++r)
              alphaE += Qaa(r, i)*abit[r];

            for (size_t ib = 0; ib < lb; ++ib) {
              const bitset<nbit__> bbit = block->string_bits_b(ib);

              double betaE = 0.0;
              for (int r = 0; r < rnorb; ++r)
                betaE += Qbb(r, i)*bbit[r];

              data_base[ib] += alphaE + betaE;
            }
          }
        }
      }
    }
  }
  denom_t.tick_print("denom: mixed");
}
