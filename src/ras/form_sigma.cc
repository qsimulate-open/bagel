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

  return operator()(ccvec, jop->mo1e()->matrix(), jop->mo2e(), conv);
}

void FormSigmaRAS::operator()(const RASCivecView cc, RASCivecView sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<const Matrix> mo2e = jop->mo2e();
  const int norb = jop->nocc();

  auto mo2e_hz = [&norb, &mo2e] (const int i, const int j, const int k, const int l) { return mo2e->element(i + norb*j, k + norb*l); };

  Matrix g(*jop->mo1e()->matrix());

  for (int k = 0, kl = 0; k < norb; ++k) {
    for (int l = 0; l < k; ++l, ++kl) {
      { // g_kl
        double val = -mo2e_hz(k, k, k, l);
          for (int j = 0; j < k; ++j) val -= mo2e_hz(k,j,j,l);
        g(l,k) += val;
      }

      { // g_lk
        double val = 0.0;
        for (int j = 0; j < l; ++j) val -= mo2e_hz(l,j,j,k);
        g(k,l) += val;
      }
    }
    // g_kk
    double val = -0.5*mo2e_hz(k,k,k,k);
    for (int j = 0; j < k; ++j) val -= mo2e_hz(k,j,j,k);
    g(k,k) += val;
    ++kl;
  }

  // (taskaa)
  sigma_aa(cc, sigma, g.data(), mo2e->data());

  // (taskbb)
  sigma_bb(cc, sigma, g.data(), mo2e->data());

  // (taskab) alpha-beta contributions
  sigma_ab(cc, sigma, mo2e->data());
}

shared_ptr<RASDvec> FormSigmaRAS::operator()(shared_ptr<const RASDvec> ccvec, shared_ptr<const Matrix> mo1e, shared_ptr<const Matrix> mo2e, const vector<int>& conv) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();
  const int norb = det->norb();

  auto mo2e_hz = [&norb, &mo2e] (const int i, const int j, const int k, const int l) { return mo2e->element(i + norb*j, k + norb*l); };

  Matrix g(norb, norb);
  if (mo1e)
    g += *mo1e;

  if (mo2e) {
    for (int k = 0, kl = 0; k < norb; ++k) {
      for (int l = 0; l < k; ++l, ++kl) {
        { // g_kl
          double val = -mo2e_hz(k, k, k, l);
            for (int j = 0; j < k; ++j) val -= mo2e_hz(k,j,j,l);
          g(l,k) += val;
        }

        { // g_lk
          double val = 0.0;
          for (int j = 0; j < l; ++j) val -= mo2e_hz(l,j,j,k);
          g(k,l) += val;
        }
      }
      // g_kk
      double val = -0.5*mo2e_hz(k,k,k,k);
      for (int j = 0; j < k; ++j) val -= mo2e_hz(k,j,j,k);
      g(k,k) += val;
      ++kl;
    }
  }

  auto sigmavec = make_shared<RASDvec>(det, nstate);

  // Not sure why you would ever do this, but just in case
  if (!mo1e && !mo2e)
    return sigmavec;

  // Bit of a temporary hack to make life easier if no mo2e is provided
  shared_ptr<const Matrix> twoelectron = ( !mo2e ? make_shared<Matrix>(norb*norb, norb*norb) : mo2e );

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
#ifdef HAVE_MPI_H
    if ( istate % mpi__->size() == mpi__->rank() ) {
#endif
      Timer pdebug(2);
      const RASCivecView cc(*ccvec->data(istate));
      RASCivecView sigma(*sigmavec->data(istate));

      // (taskaa)
      sigma_aa(cc, sigma, g.data(), twoelectron->data());
      pdebug.tick_print("taskaa");

      // (taskbb)
      sigma_bb(cc, sigma, g.data(), twoelectron->data());
      pdebug.tick_print("taskbb");

      // (taskab) alpha-beta contributions
      if (mo2e)
        sigma_ab(cc, sigma, twoelectron->data());
      pdebug.tick_print("taskab");
#ifdef HAVE_MPI_H
    }
#endif
  }

#ifdef HAVE_MPI_H
  for (int istate = 0; istate != nstate; ++istate) {
    if (!conv[istate])
      mpi__->broadcast(sigmavec->data(istate)->data(), sigmavec->data(istate)->size(), istate % mpi__->size());
  }
#endif

  return sigmavec;
}

// sigma_2 in the Olsen paper
void FormSigmaRAS::sigma_aa(const RASCivecView cc, RASCivecView sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASDeterminants> det = cc.det();
  assert(*det == *sigma.det());

  const int norb = det->norb();
  const size_t la = det->lena();

  Matrix F(la, batchsize_);

  // Let's just get it working first, thread it later
  for (auto& ispace : *det->stringspacea()) {
    // Do this multiplication batchwise
    const int nbatches = (ispace->size() - 1)/batchsize_ + 1;
    for (int batch = 0; batch < nbatches; ++batch) {
      const size_t batchstart = batch * batchsize_;
      const size_t batchlength = min(static_cast<size_t>(batchsize_), ispace->size() - batchstart);

      F.zero();

      for (size_t ia = 0; ia < batchlength; ++ia) {
        double * const fdata = F.element_ptr(0, ia);
        const size_t offset = batchstart + ispace->offset();
        for (auto& iterkl : det->phia(ia + offset)) {
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
      // S(beta, alpha) += C(beta, alpha) * F(alpha, alpha')
      for (auto& iblock : cc.blocks()) {
        if (!iblock) continue;
        if (!det->allowed(ispace, iblock->stringsb())) continue;
        shared_ptr<RASBlock<double>> target_block = sigma.block(iblock->stringsb(), ispace);

        assert(iblock->lenb() == target_block->lenb());
        assert(ispace->size() == target_block->lena());
        dgemm_("N", "N", target_block->lenb(), batchlength, iblock->lena(), 1.0,
                         iblock->data(), iblock->lenb(),
                         F.element_ptr(iblock->stringsa()->offset(), 0), F.ndim(), 1.0,
                         target_block->data() + batchstart * target_block->lenb(), target_block->lenb());
      }
    }
  }
}

void FormSigmaRAS::sigma_bb(const RASCivecView cc, RASCivecView sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASCivec> cc_trans = cc.transpose();
  auto sig_trans = make_shared<RASCivec>(cc_trans->det());

  sigma_aa(RASCivecView(*cc_trans), RASCivecView(*sig_trans), g, mo2e);

  sigma.ax_plus_y(1.0, *sig_trans->transpose(sigma.det()));
}

// Helper class to deal with C' matrices
namespace bagel {
  class Cprime {
    protected:
      // matrices named for which space (RASI, II, or III) runs first.
      // remaining spaces run in original order (II, I, III)
      shared_ptr<const RASString> source_space_;
      shared_ptr<const RASDeterminants> det_;

      // store all the information relating to which subspaces of C' are present in the matrices
      vector<pair<const DetMapBlock*, shared_ptr<Matrix>>> blocks_;

    public:
      Cprime(shared_ptr<const RASString> space, shared_ptr<const RASDeterminants> det,
          vector<pair<const DetMapBlock*, shared_ptr<Matrix>>>&& data) : source_space_(space), det_(det), blocks_(move(data)) { }

      shared_ptr<Matrix> get_matrix(shared_ptr<const RASString> target_space) const {
        const size_t stringsize = source_space_->size();
        vector<vector<size_t>> indices;
        size_t nallowed = 0;
        for (auto& b : blocks_) {
          size_t index = 0;
          vector<size_t> tmp_indices;
          for(auto& i : *b.first) {
            if (det_->allowed(target_space->strings(0), det_->string_bits_b(i.target))) {
              tmp_indices.push_back(index);
            }
            ++index;
          }
          nallowed += tmp_indices.size();
          indices.emplace_back(move(tmp_indices));
        }

        if (nallowed == 0)
        { return nullptr; }
        else {
          size_t col = 0;
          auto out = make_shared<Matrix>(stringsize, nallowed);
          auto biter = blocks_.begin();
          for (auto i = indices.begin(); i != indices.end(); ++i, ++biter) {
            if (i->size() > 0) {
              shared_ptr<Matrix> mat = biter->second;
              for (auto j : *i)
                copy_n(mat->element_ptr(0, j), stringsize, out->element_ptr(0, col++));
            }
          }
          return out;
        }
      }
  };
}


void FormSigmaRAS::sigma_ab(const RASCivecView cc, RASCivecView sigma, const double* mo2e) const {
  assert(*cc.det() == *sigma.det());
  shared_ptr<const RASDeterminants> det = cc.det();

  const int norb = det->norb();

  map<size_t, map<size_t, pair<vector<tuple<size_t, int, int>>, shared_ptr<SparseMatrix>>>> Fmatrices;

  for (auto& ispace : *det->stringspacea()) {
    const int nspaces = det->stringspacea()->nspaces();
    const size_t la = ispace->size();

    // These are for building the initial versions of the sparse matrices
    vector<vector<double>> data(nspaces);
    vector<vector<int>> cols(nspaces);
    vector<vector<int>> rind(nspaces);
    vector<vector<tuple<size_t, int, int>>> sparse_info(nspaces);

    vector<pair<size_t, int>> bounds;
    for (auto& isp : *det->stringspacea()) {
      bounds.emplace_back(isp->offset(), isp->offset() + isp->size());
    }
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
      else Fmap.emplace(bounds[isp].first, make_pair(move(sparse_info[isp]), nullptr));
    }

    Fmatrices.emplace(ispace->offset(), move(Fmap));
  }

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      // L(I), R(I), sign(I) building
      const size_t phisize = accumulate(det->phib_ij(ij).begin(), det->phib_ij(ij).end(), 0ull, [] (size_t i, const DetMapBlock& m) { return i + m.size(); });
      if (phisize == 0) continue;

      map<shared_ptr<const RASString>, shared_ptr<Cprime>> Cp_map;

      // gathering
      {
        map<shared_ptr<const RASString>, vector<pair<const DetMapBlock*, shared_ptr<Matrix>>>> Cp_temp;

        for ( auto& iphiblock : det->phib_ij(ij) ) {
          vector<shared_ptr<const RASBlock<double>>> blks = cc.allowed_blocks<1>(iphiblock.space());
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

            Cp_temp[iblock->stringsa()].emplace_back(&iphiblock, tmp);
          }
        }

        for (auto& cpblock : Cp_temp)
          Cp_map.emplace(cpblock.first, make_shared<Cprime>(cpblock.first, det, move(cpblock.second)));
      }

      // build V(I), block by block
      for (auto& ispace : *det->stringspacea()) {
        const size_t la = ispace->size();

        // build reduced version of phiblock and Cp
        vector<DetMapBlock> reduced_phi;
        size_t offset = 0;
        for (auto& phiblock : det->phib_ij(ij)) {
          vector<DetMap> phis;
          for (auto& phi : phiblock) {
            shared_ptr<const RASString> betaspace = det->space<1>(det->string_bits_b(phi.target));
            if (det->allowed(ispace, betaspace))
              phis.emplace_back(phi);
          }
          if (phis.size() > 0) {
            const size_t sz = phis.size();
            reduced_phi.emplace_back(offset, phiblock.space(), move(phis));
            offset += sz;
          }
        }
        const size_t reduced_phisize = offset;

        if (reduced_phisize == 0) continue;

        auto Vt = make_shared<Matrix>(la, reduced_phisize);

        // Grabbing sparse matrices
        auto& Fmap = Fmatrices.at(ispace->offset());
        const double* mo2e_ij = mo2e + i + norb*norb*j;

        for (auto& cpblock : Cp_map) {
          shared_ptr<const RASString> source_space = cpblock.first;
          shared_ptr<Cprime> cp = cpblock.second;
          shared_ptr<Matrix> cp_matrix = cp->get_matrix(ispace);
          if (cp_matrix) {
            shared_ptr<SparseMatrix> Ft_block = Fmap[source_space->offset()].second;
            if (Ft_block) {
              // Update data
              Ft_block->zero();
              double* fdata = Ft_block->data();
              vector<tuple<size_t, int, int>>& replace_data = Fmap[source_space->offset()].first;
              for_each(replace_data.begin(), replace_data.end(),
                        [&mo2e_ij, &fdata] (tuple<size_t, int, int> i)
                          { fdata[get<0>(i)] += static_cast<double>(get<1>(i)) * mo2e_ij[get<2>(i)]; });

              // multiply
              auto Vt_block = make_shared<Matrix>(*Ft_block * *cp_matrix);
              size_t current = 0;
              for (auto& phiblock : reduced_phi) {
                if (!det->allowed(phiblock.space(), source_space)) continue;
                Vt->add_block(1.0, 0, phiblock.offset(), Vt_block->ndim(), phiblock.size(), Vt_block->element_ptr(0, current));
                current += phiblock.size();
              }
            }
          }
        }

        // scatter
        double* vdata = Vt->data();
        for (auto& iphiblock : reduced_phi ) {
          for (auto& iphi : iphiblock) {
            shared_ptr<const RASString> betaspace = det->space<1>(det->string_bits_b(iphi.target));
            const double* sourcedata = vdata;

            shared_ptr<RASBlock<double>> sgblock = sigma.block(betaspace, ispace);
            double* targetdata = sgblock->data() + iphi.target - betaspace->offset();

            const size_t lb = sgblock->lenb();

            for (size_t i = 0; i < la; ++i, targetdata+=lb, ++sourcedata) {
              *targetdata += *sourcedata;
            }

            vdata += la;
          }
        }
      }
    }
  }
}
