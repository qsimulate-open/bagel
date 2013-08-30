//
// BAGEL - Parallel electron correlation program.
// Filename: ras/rasci_compute.cc
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

#include <src/ras/rasci.h>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <src/util/constants.h>
#include <src/util/taskqueue.h>
#include <vector>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Olsen */
vector<shared_ptr<RASCivec>> RASCI::form_sigma(const vector<shared_ptr<RASCivec>>& ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const {
  const int ij = norb_*norb_;

  const int nstate = nstate_;

  vector<shared_ptr<RASCivec>> sigmavec;
  for (int i = 0; i < nstate; ++i) { sigmavec.push_back( make_shared<RASCivec>(det_) ); }

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(0);
    if (conv[istate]) continue;
    shared_ptr<const RASCivec> cc = ccvec.at(istate);
    shared_ptr<RASCivec> sigma = sigmavec.at(istate);

    // (taskaa)
    sigma_aa(cc, sigma, jop);
    pdebug.tick_print("taskaa");

    // (taskbb)
    sigma_bb(cc, sigma, jop);
    pdebug.tick_print("taskbb");

    // (taskab) alpha-beta contributions
    sigma_ab(cc, sigma, jop);
    pdebug.tick_print("taskab");
  }

  return sigmavec;
}

// sigma_2 in the Olsen paper
void RASCI::sigma_aa(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  assert(det == sigma->det());

  const int la = det->lena();

  const int norb = norb_;
  unique_ptr<double[]> g(new double[norb*norb]);
  for (int k = 0, kl = 0; k < norb; ++k) {
    for (int l = 0; l < k; ++l, ++kl) {
      { // g_kl
        double val = jop->mo1e(kl) - jop->mo2e_hz(k, k, k, l);
        for (int j = 0; j < k; ++j) val -= jop->mo2e_hz(k,j,j,l);
        g[l + k * norb] = val;
      }

      { // g_lk
        double val = jop->mo1e(kl);
        for (int j = 0; j < l; ++j) val -= jop->mo2e_hz(l,j,j,k);
        g[k + l * norb] = val;
      }
    }
    // g_kk
    double val = jop->mo1e(kl) - 0.5*jop->mo2e_hz(k,k,k,k);
    for (int j = 0; j < k; ++j) val -= jop->mo2e_hz(k,j,j,k);
    g[k + k * norb] = val;
    ++kl;
  }

  // Let's just get it working first, thread it later
  int aiter = 0;
  for (auto& ispace : det->stringspacea()) {
    if (!ispace) continue;
    unique_ptr<double[]> F(new double[la * ispace->size()]);
    fill_n(F.get(), la * ispace->size(), 0.0);
    double* data = F.get();
    for (int ia = 0; ia < ispace->size(); ++ia, ++aiter, data+=la) {
      for (auto& iterkl : det->phia(aiter)) {
        data[iterkl.source] += static_cast<double>(iterkl.sign) * g[iterkl.ij];
        for (auto& iterij : det->phia(iterkl.source)) {
          if (iterij.ij < iterkl.ij) continue;
          const int ii = iterij.ij/norb;
          const int jj = iterij.ij%norb;
          const int kk = iterkl.ij/norb;
          const int ll = iterkl.ij%norb;
          data[iterij.source] += static_cast<double>(iterkl.sign*iterij.sign) * (iterkl.ij == iterij.ij ? 0.5 : 1.0) * jop->mo2e_hz(ii, kk, jj, ll);
        }
      }
    }

    // F is finished, matrix-matrix multiply (but to the right place)
    for (auto& iblock : cc->blocks()) {
      if (!iblock) continue;
      shared_ptr<RASBlock<double>> target_block = sigma->block(ispace->nholes(), iblock->stringb()->nholes(), ispace->nparticles(), iblock->stringb()->nparticles());
      if (!target_block) continue;

      dgemm_("N", "N", iblock->lenb(), ispace->size(), iblock->lena(), 1.0, iblock->data(), iblock->lenb(),
        F.get() + iblock->stringa()->offset(), la, 1.0, target_block->data(), target_block->lenb());
    }
  }
}

void RASCI::sigma_bb(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<const RASCivec> cc_trans = cc->transpose();
  auto sig_trans = make_shared<RASCivec>(cc_trans->det());

  sigma_aa(cc_trans, sig_trans, jop);

  sigma->ax_plus_y(1.0, *sig_trans->transpose(sigma->det()));
}

void RASCI::sigma_ab(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, shared_ptr<const MOFile> jop) const {
  assert(cc->det() == sigma->det());
  shared_ptr<const RASDeterminants> det = cc->det();

  const int norb = norb_;

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      vector<tuple<bitset<nbit__>, bitset<nbit__>, int>> philist;
      // L(I), R(I), sign(I) building
      for (auto& bbit : det->stringb()) {
        if ( !bbit[i] || ( i!=j && bbit[j]) ) continue;
        bitset<nbit__> sourcebbit = bbit; sourcebbit.reset(i); sourcebbit.set(j);
        if ( !det->allowed(sourcebbit) ) continue;
        int sign = det->sign(bbit, i, j);
        philist.emplace_back(sourcebbit, bbit, sign);
      }
      if (philist.size() == 0) continue;

      unique_ptr<double[]> Cp(new double[philist.size() * det->lena()]);
      fill_n(Cp.get(), philist.size() * det->lena(), 0.0);
      double* cdata = Cp.get();

      // gathering
      for ( auto& abit : det->stringa() ) {
        for ( auto& iphi : philist ) {
          // Cp(I, Ja) = sign(I) * C(L(I), Ja)
          if ( det->allowed(abit, get<0>(iphi)) ) *cdata = cc->element(get<0>(iphi), abit) * static_cast<double>(get<2>(iphi));
          ++cdata;
        }
      }

      // build F, block by block
      for (auto& ispace : det->stringspacea()) {
        if (!ispace) continue;
        const int la = ispace->size();

        unique_ptr<double[]> F(new double[ la * det->lena() ]);
        double* fdata = F.get();
        fill_n(fdata, la * det->lena(), 0.0);
        for (int ia = 0; ia < la; ++ia, fdata+=det->lena()) {
          for (auto& iter : det->phia(ia + ispace->offset())) {
            const int kk = iter.ij/norb;
            const int ll = iter.ij%norb;
            fdata[iter.source] += static_cast<double>(iter.sign) * jop->mo2e_hz(i, kk, j, ll);
          }
        }

        unique_ptr<double[]> VI(new double[ la * philist.size() ]);
        dgemm_("N", "N", philist.size(), la, det->lena(), 1.0, Cp.get(), philist.size(), F.get(), det->lena(), 0.0, VI.get(), philist.size());

        // scatter
        double* vdata = VI.get();
        for ( auto& abit : *ispace ) {
          for (auto& iphi : philist) {
            if (det->allowed(abit, get<1>(iphi))) sigma->element(get<1>(iphi), abit) += *vdata;
            ++vdata;
          }
        }
      }
    }
  }
}
