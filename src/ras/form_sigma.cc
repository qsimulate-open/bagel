//
// BAGEL - Parallel electron correlation program.
// Filename: ras/form_sigma.cc
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

#include <src/ras/form_sigma.h>
#include <src/math/sparsematrix.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Olsen */
// This is how this is accessed from RASCI
shared_ptr<RASDvec> FormSigmaRAS::operator()(shared_ptr<const RASDvec> ccvec, shared_ptr<const MOFile> jop,
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

  auto sigmavec = make_shared<RASDvec>(det, nstate);

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(2);
    if (conv[istate]) continue;
    shared_ptr<const RASCivec> cc = ccvec->data(istate);
    shared_ptr<RASCivec> sigma = sigmavec->data(istate);

    // (taskaa)
    sigma_aa(cc, sigma, g.get(), jop->mo2e_ptr());
    pdebug.tick_print("taskaa");

    // (taskbb)
    sigma_bb(cc, sigma, g.get(), jop->mo2e_ptr());
    pdebug.tick_print("taskbb");

    // (taskab) alpha-beta contributions
    sigma_ab(cc, sigma, jop->mo2e_ptr());
    pdebug.tick_print("taskab");
  }

  return sigmavec;
}

// This is how this is accessed from MEH
shared_ptr<RASDvec> FormSigmaRAS::operator()(shared_ptr<const RASDvec> ccvec, const double* mo1e, const double* mo2e) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();
  const int norb = det->norb();

  auto mo2e_hz = [&mo2e, &norb] (const int i, const int j, const int k, const int l) { return mo2e[i + j*norb + norb*norb*(k + l*norb)]; };

  unique_ptr<double[]> g(new double[norb*norb]);
  for (int k = 0, kl = 0; k < norb; ++k) {
    for (int l = 0; l < k; ++l, ++kl) {
      { // g_kl
        double val = mo1e[l + k * norb] - mo2e_hz(k, k, k, l);
        for (int j = 0; j < k; ++j) val -= mo2e_hz(k,j,j,l);
        g[l + k * norb] = val;
      }

      { // g_lk
        double val = mo1e[k + l * norb];
        for (int j = 0; j < l; ++j) val -= mo2e_hz(l,j,j,k);
        g[k + l * norb] = val;
      }
    }
    // g_kk
    double val = mo1e[k + k*norb] - 0.5*mo2e_hz(k,k,k,k);
    for (int j = 0; j < k; ++j) val -= mo2e_hz(k,j,j,k);
    g[k + k * norb] = val;
    ++kl;
  }

  auto sigmavec = make_shared<RASDvec>(det, nstate);

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(0);
    shared_ptr<const RASCivec> cc = ccvec->data(istate);
    shared_ptr<RASCivec> sigma = sigmavec->data(istate);

    // (taskaa)
    sigma_aa(cc, sigma, g.get(), mo2e);
    pdebug.tick_print("taskaa");

    // (taskbb)
    sigma_bb(cc, sigma, g.get(), mo2e);
    pdebug.tick_print("taskbb");

    // (taskab) alpha-beta contributions
    sigma_ab(cc, sigma, mo2e);
    pdebug.tick_print("taskab");
  }

  return sigmavec;
}

// A bit of a temporary hack for 1e terms
shared_ptr<RASDvec> FormSigmaRAS::operator()(shared_ptr<const RASDvec> ccvec, const double* mo1e) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();
  const int norb = det->norb();

  auto sigmavec = make_shared<RASDvec>(det, nstate);

  unique_ptr<double[]> blank2e(new double[norb*norb*norb*norb]);
  fill_n(blank2e.get(), norb*norb*norb*norb, 0.0);

  for (int istate = 0; istate != nstate; ++istate) {
    shared_ptr<const RASCivec> cc = ccvec->data(istate);
    shared_ptr<RASCivec> sigma = sigmavec->data(istate);

    // (taskaa)
    sigma_aa(cc, sigma, mo1e, blank2e.get());

    // (taskbb)
    sigma_bb(cc, sigma, mo1e, blank2e.get());
  }

  return sigmavec;
}

// sigma_2 in the Olsen paper
void FormSigmaRAS::sigma_aa(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  assert(*det == *sigma->det());

  const int norb = det->norb();
  const size_t la = det->lena();

  // Let's just get it working first, thread it later
  for (auto& ispace : det->stringspacea()) {
    if (!ispace) continue;
    unique_ptr<double[]> F(new double[la * ispace->size()]);
    fill_n(F.get(), la * ispace->size(), 0.0);
    double* fdata = F.get();
    for (size_t ia = 0; ia < ispace->size(); ++ia, fdata+=la) {
      for (auto& iterkl : det->phia(ia + ispace->offset())) {
        fdata[iterkl.source] += static_cast<double>(iterkl.sign) * g[iterkl.ij];
        for (auto& iterij : det->phia(iterkl.source)) {
          if (iterij.ij < iterkl.ij) continue;
          const int ii = iterij.ij/norb;
          const int jj = iterij.ij%norb;
          const int kk = iterkl.ij/norb;
          const int ll = iterkl.ij%norb;
          fdata[iterij.source] += static_cast<double>(iterkl.sign*iterij.sign) * (iterkl.ij == iterij.ij ? 0.5 : 1.0) * mo2e[ii + kk*norb + norb*norb*(jj + ll * norb)];
        }
      }
    }

    // F is finished, matrix-matrix multiply (but to the right place)
    for (auto& iblock : cc->blocks()) {
      if (!iblock) continue;
      if (!det->allowed(ispace, iblock->stringb())) continue;
      shared_ptr<RASBlock<double>> target_block = sigma->block(iblock->stringb(), ispace);

      assert(iblock->lenb() == target_block->lenb());
      assert(ispace->size() == target_block->lena());
      dgemm_("N", "N", target_block->lenb(), target_block->lena(), iblock->lena(), 1.0, iblock->data(), iblock->lenb(),
        F.get() + iblock->stringa()->offset(), la, 1.0, target_block->data(), target_block->lenb());
    }
  }
}

void FormSigmaRAS::sigma_bb(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASCivec> cc_trans = cc->transpose();
  auto sig_trans = make_shared<RASCivec>(cc_trans->det());

  sigma_aa(cc_trans, sig_trans, g, mo2e);

  sigma->ax_plus_y(1.0, *sig_trans->transpose(sigma->det()));
}

void FormSigmaRAS::sigma_ab(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* mo2e) const {
  assert(cc->det() == sigma->det());
  shared_ptr<const RASDeterminants> det = cc->det();

  const int norb = det->norb();

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

      // build F, block by block
      for (auto& ispace : det->stringspacea()) {
        if (!ispace) continue;
        const size_t la = ispace->size();

        auto Vt = make_shared<Matrix>(la, phisize);
        if ( sparse_ ) {
          // Making several SparseMatrix objects
          // First, how many?
          const int nspaces = accumulate(det->stringspacea().begin(), det->stringspacea().end(), 0, [] (int i, shared_ptr<const StringSpace> s) { return i + ( s ? 1 : 0); });
          vector<vector<double>> data(nspaces, vector<double>());
          vector<vector<int>> cols(nspaces, vector<int>());
          vector<vector<int>> rind(nspaces, vector<int>());

          vector<pair<size_t, int>> bounds;
          for (auto& isp : det->stringspacea()) if (isp) bounds.emplace_back(isp->offset(), isp->offset() + isp->size());
          assert(bounds.size() == nspaces);

          // figures out which space a column belongs to and subtracts the offset
          auto which_space = [&bounds] (int& col) {
            int space = 0;
            for (auto& b : bounds) {
              if ( col >= b.first && col < b.second ) {
                col -= b.first;
                return space;
              }
              ++space;
            }
            assert(false); return 0;
          };

          // Not sure how to guess how many elements would be filled in each stringspace

          for (int ia = 0; ia < la; ++ia) {
            vector<map<size_t, double>> rows(nspaces);
            for (auto& iter : det->phia(ia + ispace->offset())) {
              const int kk = iter.ij/norb;
              const int ll = iter.ij%norb;
              int col = iter.source;
              const int sp = which_space(col);
              rows.at(sp)[col] += static_cast<double>(iter.sign) * mo2e[i + norb*kk + norb*norb*(j + norb*ll)];
            }

            for (int isp = 0; isp < nspaces; ++isp) {
              rind.at(isp).push_back(data[isp].size() + 1);
              for (auto& irow : rows.at(isp)) {
                cols[isp].push_back(irow.first + 1);
                data[isp].push_back(irow.second);
              }
            }
          }

          map<size_t, shared_ptr<SparseMatrix>> Ft_map;

          for (int isp = 0; isp < nspaces; ++isp) {
            if (data.at(isp).size() > 0) {
              rind.at(isp).push_back(data.at(isp).size() + 1);
              const int mdim = bounds.at(isp).second - bounds.at(isp).first;
              Ft_map.emplace(bounds.at(isp).first, make_shared<SparseMatrix>(la, mdim, data.at(isp), cols.at(isp), rind.at(isp)));
            }
            else Ft_map.emplace(bounds.at(isp).first, shared_ptr<SparseMatrix>());
          }

          for (auto& iphiblock : det->phib_ij(ij)) {
            vector<shared_ptr<const StringSpace>> allowed_spaces = det->allowed_spaces<1>(iphiblock.space());
            for (auto& mult_space : allowed_spaces) {
              shared_ptr<Matrix> Cp_block = Cp_map.at(make_pair(mult_space->offset(), iphiblock.offset()));
              shared_ptr<SparseMatrix> Ft_block = Ft_map.at(mult_space->offset());
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
}
