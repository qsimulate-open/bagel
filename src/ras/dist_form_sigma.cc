//
// BAGEL - Parallel electron correlation program.
// Filename: ras/dist_form_sigma.cc
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

#include <map>

#include <src/ras/dist_form_sigma.h>
#include <src/math/sparsematrix.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Olsen */
// This is how this is accessed from RASCI
shared_ptr<DistRASDvec> DistFormSigmaRAS::operator()(shared_ptr<const DistRASDvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();

  const int norb = det->norb();
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

  auto sigmavec = make_shared<DistRASDvec>(det, nstate);

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(2);
    if (conv[istate]) continue;
    shared_ptr<const DistRASCivec> cc = ccvec->data(istate);
    shared_ptr<DistRASCivec> sigma = sigmavec->data(istate);

    // start transpose
    shared_ptr<const DistRASCivec> cctrans = cc->transpose();

    // (taskbb)
    sigma_bb(cc, sigma, g.get(), jop->mo2e_ptr());
    pdebug.tick_print("taskbb");

    // finish transpose
    cctrans->transpose_wait();
    pdebug.tick_print("wait");

    // (taskaa)
    shared_ptr<DistRASCivec> strans = cctrans->clone();
    sigma_bb(cctrans, strans, g.get(), jop->mo2e_ptr());
    pdebug.tick_print("taskaa");

    // start transpose back
    shared_ptr<const DistRASCivec> saa = strans->transpose();

    // (taskab) alpha-beta contributions
    sigma_ab(cc, sigma, jop->mo2e_ptr());
    pdebug.tick_print("taskab");

    // finish transpose back
    saa->transpose_wait();
    pdebug.tick_print("wait1");
    sigma->ax_plus_y(1.0, *saa);
  }

  return sigmavec;
}

// A bit of a temporary hack for 1e terms
shared_ptr<DistRASDvec> DistFormSigmaRAS::operator()(shared_ptr<const DistRASDvec> ccvec, const double* mo1e) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();
  const int norb = det->norb();

  auto sigmavec = make_shared<DistRASDvec>(det, nstate);

  unique_ptr<double[]> blank2e(new double[norb*norb*norb*norb]);
  fill_n(blank2e.get(), norb*norb*norb*norb, 0.0);

  for (int istate = 0; istate != nstate; ++istate) {
    shared_ptr<const DistRASCivec> cc = ccvec->data(istate);
    shared_ptr<DistRASCivec> sigma = sigmavec->data(istate);

    // start transpose
    shared_ptr<const DistRASCivec> cctrans = cc->transpose();

    // (taskbb)
    sigma_bb(cc, sigma, mo1e, blank2e.get());

    shared_ptr<DistRASCivec> strans = cctrans->clone();
    cctrans->transpose_wait();

    // (taskaa)
    sigma_bb(cctrans, strans, mo1e, blank2e.get());

    shared_ptr<const DistRASCivec> saa = strans->transpose();
    saa->transpose_wait();

    sigma->ax_plus_y(1.0, *saa);
  }

  return sigmavec;
}

// sigma_bb requires no communication
void FormSigmaRAS::sigma_bb(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  assert(*det = *sigma->det());

  const int norb = det->norb();
  const size_t lb = det->lenb();

  for (auto& ispace : det->stringspaceb()) {
    if (!ispace) continue;
    auto F = make_shared<Matrix>(lb, ispace->size());
    double* fdata = F->data();
    for (size_t ib = 0; ib < ispace->size(); ++ib, fdata+=lb) {
      for (auto& iterkl : det->phib(ib + ispace->offset())) {
        fdata[iterkl.source] += static_cast<double>(iterkl.sign) * g[iterkl.ij];
        for (auto& iterij : det->phib(iterkl.source)) {
          if (iterij.ij < iterkl.ij) continue
          const int ii = iterij.ij/norb;
          const int jj = iterij.ij%norb;
          const int kk = iterkl.ij/norb;
          const int ll = iterkl.ij%norb;
          fdata[iterij.source] += static_cast<double>(iterkl.sign * iterij.sign) * (iterkl.ij == iterij.ij ? 0.5 : 1.0) * mo2e[ii + kk*norb + norb*norb*(jj + ll * norb)];
        }
      }
    }

    // F is finished, apply to the right places
    for (auto& iblock : cc->blocks()) {
      if (!iblock) continue;
      if (!det->allowed(ispace, iblock->stringb())) continue;
      shared_ptr<DistRASBlock<double>> tblock = sigma->block(iblock->stringb(), ispace);

      dgemm_("T", "N", ispace->size(), iblock->asize(), iblock->lenb(), 1.0,
        F->element_ptr(0, iblock->stringb()->offset()), F->ndim(), iblock->local(), iblock->lenb(), 1.0, tblock->local(), tblock->lenb());
    }
  }
}

void FormSigmaRAS::sigma_ab(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* mo2e) const {
#if 0
  assert(cc->det() == sigma->det());
  shared_ptr<const RASDeterminants> det = cc->det();

  const int norb = det->norb();

  map<size_t, map<size_t, pair<vector<tuple<size_t, int, int>>, shared_ptr<SparseMatrix>>>> Fmatrices;

  if (sparse_) {
    for (auto& ispace : det->stringspacea()) {
      if (!ispace) continue;
      const int nspaces = accumulate(det->stringspacea().begin(), det->stringspacea().end(), 0, [] (int i, shared_ptr<const StringSpace> s) { return i + ( s ? 1 : 0); });
      const size_t la = ispace->size();

      // These are for building the initial versions of the sparse matrices
      vector<vector<double>> data(nspaces);
      vector<vector<int>> cols(nspaces);
      vector<vector<int>> rind(nspaces);
      vector<vector<tuple<size_t, int, int>>> sparse_info(nspaces);

      vector<pair<size_t, int>> bounds;
      for (auto& isp : det->stringspacea()) if (isp) bounds.emplace_back(isp->offset(), isp->offset() + isp->size());
      assert(bounds.size() == nspaces);

      for (int ia = 0; ia < la; ++ia) {
        map<size_t, double> row;
        map<size_t, vector<tuple<int, int>>> row_positions;
        for (auto& iter : det->phia(ia + ispace->offset())) {
          const int kk = iter.ij/norb;
          const int ll = iter.ij%norb;
          row_positions[iter.source].emplace_back(iter.sign, ll*norb + kk*norb*norb*norb);
        }

        for (int sp = 0; sp < nspaces; ++sp) rind[sp].push_back(data[sp].size() + 1);

        auto ibound = bounds.begin();
        int sp = 0;
        for (auto& irow : row_positions) {
          while (irow.first >= ibound->second) { ++ibound; ++sp; }
          const int pos = data[sp].size();
          for (auto& i : irow.second)
            sparse_info[sp].emplace_back(pos, get<0>(i), get<1>(i));
          cols[sp].push_back(irow.first + 1 - ibound->first);
          data[sp].push_back(1.0);
        }
      }

      map<size_t, pair<vector<tuple<size_t, int, int>>, shared_ptr<SparseMatrix>>> Fmap;

      for (int isp = 0; isp < nspaces; ++isp) {
        if (data[isp].size() > 0) {
          rind[isp].push_back(data[isp].size() + 1);
          const int mdim = bounds.at(isp).second - bounds.at(isp).first;
          Fmap.emplace(bounds[isp].first, make_pair(move(sparse_info[isp]), make_shared<SparseMatrix>(la, mdim, data[isp], cols[isp], rind[isp])));
        }
        else Fmap.emplace(bounds[isp].first, make_pair(move(sparse_info[isp]), shared_ptr<SparseMatrix>()));
      }

      Fmatrices.emplace(ispace->offset(), move(Fmap));
    }
  }

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      // L(I), R(I), sign(I) building
      const size_t phisize = accumulate(det->phib_ij(ij).begin(), det->phib_ij(ij).end(), 0ull, [] (size_t i, const RAS::DMapBlock& m) { return i + m.size(); });
      if (phisize == 0) continue;

      map<pair<size_t, size_t>, shared_ptr<Matrix>> Cp_map;

      // gathering
      for ( auto& iphiblock : det->phib_ij(ij) ) {
        vector<shared_ptr<const RASBlock<double>>> blks = cc->allowed_blocks<1>(iphiblock.space());
        for (auto& iblock : blks) {
          auto tmp = make_shared<Matrix>(iblock->lena(), iphiblock.size());
          double* targetdata = tmp->data();

          const size_t lb = iblock->lenb();

          for (auto& iphi : iphiblock) {
            double sign = static_cast<double>(iphi.sign);
            const double* sourcedata = iblock->data() + iphi.source;

            for (size_t i = 0; i < iblock->lena(); ++i, ++targetdata, sourcedata+=lb)
              *targetdata = *sourcedata * sign;
          }

          Cp_map.emplace(make_pair(iblock->stringa()->offset(), iphiblock.offset()), tmp);
        }
      }

      // build V(I), block by block
      for (auto& ispace : det->stringspacea()) {
        if (!ispace) continue;
        const size_t la = ispace->size();

        auto Vt = make_shared<Matrix>(la, phisize);
        if ( sparse_ ) {
          // Making several SparseMatrix objects
          // First, how many?
          const int nspaces = accumulate(det->stringspacea().begin(), det->stringspacea().end(), 0, [] (int i, shared_ptr<const StringSpace> s) { return i + ( s ? 1 : 0); });
          vector<vector<double>> data(nspaces);
          vector<vector<int>> cols(nspaces);
          vector<vector<int>> rind(nspaces);

          vector<pair<size_t, int>> bounds;
          for (auto& isp : det->stringspacea()) if (isp) bounds.emplace_back(isp->offset(), isp->offset() + isp->size());
          assert(bounds.size() == nspaces);

          auto& Fmap = Fmatrices.at(ispace->offset());

          // Replace data in SparseMatrix
          const double* mo2e_ij = mo2e + i + norb*norb*j;
          for (auto& f : Fmap) {
            shared_ptr<SparseMatrix> sparse = f.second.second;
            if (sparse) {
              sparse->zero();
              double* fdata = sparse->data();
              for (auto& i : f.second.first)
                fdata[get<0>(i)] += static_cast<double>(get<1>(i)) * mo2e_ij[get<2>(i)];
            }
          }

          for (auto& iphiblock : det->phib_ij(ij)) {
            vector<shared_ptr<const StringSpace>> allowed_spaces = det->allowed_spaces<1>(iphiblock.space());
            for (auto& mult_space : allowed_spaces) {
              shared_ptr<Matrix> Cp_block = Cp_map.at(make_pair(mult_space->offset(), iphiblock.offset()));
              shared_ptr<SparseMatrix> Ft_block = Fmap[mult_space->offset()].second;
              if (Ft_block) {
                auto Vt_block = make_shared<Matrix>(*Ft_block * *Cp_block);
                Vt->add_block(1.0, 0, iphiblock.offset(), Vt_block->ndim(), Vt_block->mdim(), Vt_block);
              }
            }
          }
        }
        else { // dense matrix multiply
          auto F = make_shared<Matrix>( det->lena(), la );
          double* fdata = F->data();
          for (size_t ia = 0; ia < la; ++ia, fdata+=det->lena()) {
            for (auto& iter : det->phia(ia + ispace->offset())) {
              const int kk = iter.ij/norb;
              const int ll = iter.ij%norb;
              fdata[iter.source] += static_cast<double>(iter.sign) * mo2e[i + norb*kk + norb*norb*(j + norb*ll)];
            }
          }
          for (auto& iphiblock : det->phib_ij(ij)) {
            vector<shared_ptr<const StringSpace>> allowed_spaces = det->allowed_spaces<1>(iphiblock.space());
            for (auto& mult_space : allowed_spaces) {
              shared_ptr<Matrix> Cp_block = Cp_map.at(make_pair(mult_space->offset(), iphiblock.offset()));
              assert(mult_space->size() == Cp_block->ndim());
              dgemm_("T", "N", la, iphiblock.size(), mult_space->size(), 1.0, F->element_ptr(mult_space->offset(), 0), det->lena(),
                                   Cp_block->data(), Cp_block->ndim(), 1.0, Vt->element_ptr(0, iphiblock.offset()), la);
            }
          }
        }

        // scatter
        double* vdata = Vt->data();
        for (auto& iphiblock : det->phib_ij(ij) ) {
          for (auto& iphi : iphiblock) {
            shared_ptr<const StringSpace> betaspace = det->space<1>(det->stringb(iphi.target));
            if (det->allowed(ispace, betaspace)) {
              const double* sourcedata = vdata;

              shared_ptr<RASBlock<double>> sgblock = sigma->block(betaspace, ispace);
              double* targetdata = sgblock->data() + iphi.target - betaspace->offset();

              const size_t lb = sgblock->lenb();

              for (size_t i = 0; i < la; ++i, targetdata+=lb, ++sourcedata) {
                *targetdata += *sourcedata;
              }
            }

            vdata += la;
          }
        }
      }
    }
  }
#endif
}
