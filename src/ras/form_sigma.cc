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
  Matrix g(norb, norb);
  for (int k = 0, kl = 0; k < norb; ++k) {
    for (int l = 0; l < k; ++l, ++kl) {
      { // g_kl
        double val = jop->mo1e(kl) - jop->mo2e_hz(k, k, k, l);
        for (int j = 0; j < k; ++j) val -= jop->mo2e_hz(k,j,j,l);
        g(l,k) = val;
      }

      { // g_lk
        double val = jop->mo1e(kl);
        for (int j = 0; j < l; ++j) val -= jop->mo2e_hz(l,j,j,k);
        g(k,l) = val;
      }
    }
    // g_kk
    double val = jop->mo1e(kl) - 0.5*jop->mo2e_hz(k,k,k,k);
    for (int j = 0; j < k; ++j) val -= jop->mo2e_hz(k,j,j,k);
    g(k,k) = val;
    ++kl;
  }

  auto sigmavec = make_shared<RASDvec>(det, nstate);

#ifdef HAVE_MPI_H
  vector<int> requests;
#endif

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
#ifdef HAVE_MPI_H
    if ( istate % mpi__->size() == mpi__->rank() ) {
#endif
      Timer pdebug(2);
      shared_ptr<const RASCivec> cc = ccvec->data(istate);
      shared_ptr<RASCivec> sigma = sigmavec->data(istate);

      // (taskaa)
      sigma_aa(cc, sigma, g.data(), jop->mo2e_ptr());
      pdebug.tick_print("taskaa");

      // (taskbb)
      sigma_bb(cc, sigma, g.data(), jop->mo2e_ptr());
      pdebug.tick_print("taskbb");

      // (taskab) alpha-beta contributions
      sigma_ab(cc, sigma, jop->mo2e_ptr());
      pdebug.tick_print("taskab");
#ifdef HAVE_MPI_H
    }
    requests.push_back(mpi__->ibroadcast(sigmavec->data(istate)->data(), sigmavec->data(istate)->size(), istate % mpi__->size()));
#endif
  }

#ifdef HAVE_MPI_H
  for (auto& r : requests)
    mpi__->wait(r);
#endif

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
  for (auto& spaceiter : det->stringspacea()) {
    shared_ptr<const StringSpace> ispace = spaceiter.second;
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

// Helper class to deal with C' matrices
namespace bagel {
  class Cprime {
    protected:
      // matrices named for which space (RASI, II, or III) runs first.
      // remaining spaces run in original order (II, I, III)
      shared_ptr<const StringSpace> source_space_;
      shared_ptr<const RASDeterminants> det_;

      // store all the information relating to which subspaces of C' are present in the matrices
      vector<pair<const RAS::DMapBlock*, shared_ptr<Matrix>>> blocks_;

    public:
      Cprime(shared_ptr<const StringSpace> space, shared_ptr<const RASDeterminants> det,
          vector<pair<const RAS::DMapBlock*, shared_ptr<Matrix>>>&& data) : source_space_(space), det_(det), blocks_(data) { }

      shared_ptr<Matrix> get_matrix(shared_ptr<const StringSpace> target_space) const {
        const size_t stringsize = source_space_->size();
        vector<vector<size_t>> indices;
        size_t nallowed = 0;
        for (auto& b : blocks_) {
          size_t index = 0;
          vector<size_t> tmp_indices;
          for(auto& i : *b.first) {
            if (det_->allowed(target_space->strings(0), det_->stringb(i.target))) {
              tmp_indices.push_back(index);
            }
            ++index;
          }
          nallowed += tmp_indices.size();
          indices.emplace_back(move(tmp_indices));
        }

        if (nallowed == 0)
        { return shared_ptr<Matrix>(); }
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


void FormSigmaRAS::sigma_ab(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* mo2e) const {
  assert(cc->det() == sigma->det());
  shared_ptr<const RASDeterminants> det = cc->det();

  const int norb = det->norb();

  map<size_t, map<size_t, pair<vector<tuple<size_t, int, int>>, shared_ptr<SparseMatrix>>>> Fmatrices;

  for (auto& spaceiter : det->stringspacea()) {
    shared_ptr<const StringSpace> ispace = spaceiter.second;
    const int nspaces = det->stringspacea().size();
    const size_t la = ispace->size();

    // These are for building the initial versions of the sparse matrices
    vector<vector<double>> data(nspaces);
    vector<vector<int>> cols(nspaces);
    vector<vector<int>> rind(nspaces);
    vector<vector<tuple<size_t, int, int>>> sparse_info(nspaces);

    vector<pair<size_t, int>> bounds;
    for (auto& spaceiter : det->stringspacea()) {
      shared_ptr<const StringSpace> isp = spaceiter.second;
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
      else Fmap.emplace(bounds[isp].first, make_pair(move(sparse_info[isp]), shared_ptr<SparseMatrix>()));
    }

    Fmatrices.emplace(ispace->offset(), move(Fmap));
  }

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      // L(I), R(I), sign(I) building
      const size_t phisize = accumulate(det->phib_ij(ij).begin(), det->phib_ij(ij).end(), 0ull, [] (size_t i, const RAS::DMapBlock& m) { return i + m.size(); });
      if (phisize == 0) continue;

      map<shared_ptr<const StringSpace>, shared_ptr<Cprime>> Cp_map;

      // gathering
      {
        map<shared_ptr<const StringSpace>, vector<pair<const RAS::DMapBlock*, shared_ptr<Matrix>>>> Cp_temp;

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

            Cp_temp[iblock->stringa()].emplace_back(&iphiblock, tmp);
          }
        }

        for (auto& cpblock : Cp_temp)
          Cp_map.emplace(cpblock.first, make_shared<Cprime>(cpblock.first, det, move(cpblock.second)));
      }

      // build V(I), block by block
      for (auto& spaceiter : det->stringspacea()) {
        shared_ptr<const StringSpace> ispace = spaceiter.second;
        const size_t la = ispace->size();

        // build reduced version of phiblock and Cp
        vector<RAS::DMapBlock> reduced_phi;
        size_t offset = 0;
        for (auto& phiblock : det->phib_ij(ij)) {
          vector<RAS::DMap> phis;
          for (auto& phi : phiblock) {
            shared_ptr<const StringSpace> betaspace = det->space<1>(det->stringb(phi.target));
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
          shared_ptr<const StringSpace> source_space = cpblock.first;
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
            shared_ptr<const StringSpace> betaspace = det->space<1>(det->stringb(iphi.target));
            const double* sourcedata = vdata;

            shared_ptr<RASBlock<double>> sgblock = sigma->block(betaspace, ispace);
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0

// Helper class to store the differently strided versions of C'
namespace bagel {
  class Cprime {
    protected:
      // matrices named for which space (RASI, II, or III) runs first.
      // remaining spaces run in original order (II, I, III)
      array<shared_ptr<Matrix>, 3> strided_;
      shared_ptr<Matrix> trans_;
      shared_ptr<const StringSpace> space_;

      array<int, 3> nblocks_;

      // store all the information relating to which subspaces of C' are present in the matrices
      vector<const DMapBlock*> dmap_;

    public:
      shared_ptr<const Matrix> get(const int i) const { return strided_[i]; }
      shared_ptr<const Matrix> trans() const { return trans_; }
      const size_t nblocks(const int i) const { return nblocks_[i]; }

      Cprime(shared_ptr<const StringSpace> space, vector<pair<DMapBlock*, shared_ptr<Matrix>>> data) : space_(space) {
        assert(!block_data.empty());
        const size_t blocksize = accumulate(data.begin(), data.end(), 0ull,
            [] (size_t i, pair<const DMapBlock*, shared_ptr<Matrix>> p) { return i + p.second->ndim(); });

        auto r2 = make_shared<Matrix>(blocksize, block_size);
        { // II runs first initially
          int k = 0;
          for (auto& b : cpblock.second) {
            const int nn = b.second->ndim()
              r2->copy_block(k, 0, nn, block_size, *b.second);
            dmap_.push_back(b.first);
            k += nn;
          }
        }
        trans_ = r2->transpose();

        strided_[1] = r2;

        const size_t nr1 = space_->size1();
        const size_t nr2 = space_->size2();
        const size_t nr3 = space_->size3();

        nblocks_ = {{ nr2*nr3, nr1*nr3, nr1*nr2 }};

        // r1 and r3 are formed by permuting columns of r2
        if (nr1 == 1) {
          strided_[0] = r2;
        }
        else {
          shared_ptr<Matrix> r1 = r2->clone();
          for(size_t i3 = 0, k = 0; i3 < nr3; ++i3) {
            for(size_t i2 = 0; i2 < nr2; ++i2) {
              for(size_t i1 = 0; i1 < nr1; ++i1, ++k) {
                const size_t index = i2 + i1*nr2 + i3*nr1*nr2;
                copy_n(r2->element_ptr(0, index), block_size, r1->element_ptr(0, k));
              }
            }
          }
          strided_[0] = r1;
        }

        if (nr1 * nr2 == 1) {
          strided_[2] = r2;
        }
        else {
          shared_ptr<Matrix> r3 = r2->clone();
          for(size_t i1 = 0, k = 0; i1 < nr1; ++i1) {
            for(size_t i2 = 0; i2 < nr2; ++i2) {
              for(size_t i3 = 0; i3 < nr3; ++i3, ++k) {
                const size_t index = i2 + i1*nr2 + i3*nr1*nr2;
                copy_n(r2->element_ptr(0, index), block_size, r3->element_ptr(0, k));
              }
            }
          }
          strided_[2] = r3;
        }
      }
  };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FormSigmaRAS::sigma_ab_1(shared_ptr<const RASCivec> cc, shared_ptr<RASCivec> sigma, const double* mo2e) const {
  assert(cc->det() == sigma->det());
  shared_ptr<const RASDeterminants> det = cc->det();

  const int norb = det->norb();

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      // L(I), R(I), sign(I) building
      const size_t phisize = accumulate(det->phib_ij(ij).begin(), det->phib_ij(ij).end(), 0ull, [] (size_t i, const RAS::DMapBlock& m) { return i + m.size(); });
      if (phisize == 0) continue;

      const double* mo2e_ij = mo2e + i + norb*norb*j;

      map<shared_ptr<const StringSpace>, shared_ptr<Cprime>> Cp_map;

      // gathering
      {
        map<shared_ptr<const StringSpace>, vector<pair<const DMapBlock*, shared_ptr<Matrix>>>> Cp_tmp;

        // form matrices for each filled block
        for ( auto& iphiblock : det->phib_ij(ij) ) {
          vector<shared_ptr<const RASBlock<double>>> blks = cc->allowed_blocks<1>(iphiblock.space());
          for (auto& iblock : blks) {
            auto tmp = make_shared<Matrix>(iphiblock.size(), iblock->lena());
            double* targetdata = tmp->data();

            const size_t lb = iblock->lenb();

            for (auto& iphi : iphiblock) {
              double sign = static_cast<double>(iphi.sign);
              const double* sourcedata = iblock->data() + iphi.source;

              for (size_t i = 0; i < iblock->lena(); ++i)
                targetdata[i*lb] = sourcedata[i*lb] * sign;
            }

            Cp_tmp[iblock->stringa()].emplace_back(&iphiblock, tmp);
          }
        }

        // combine filled blocks corresponding to a single alpha space
        for (auto& cpblock : Cp_tmp)
          Cp_map.emplace(cpblock.first, make_shared<Cprime>(cpblock.first, cpblock.second));
      }

      // build V(I), block by block
      for (auto& spaceiter : det->stringspacea()) {
        shared_ptr<const StringSpace> ispace = spaceiter.second;
        const size_t la = ispace->size();

        auto VI_out = make_shared<Matrix>(phisize, la);
        array<shared_ptr<Matrix>, 3> VI{{ VI_out->clone(), VI_out, VI_out->clone() }};
        auto VI_trans = make_shared<Matrix>(la, phisize);
        const int nspaces = det->stringspacea().size();

        // spaces are I:0, II:1, II:2
        for(auto& cpblock : Cp_map) {
          shared_ptr<const StringSpace> source_space = cpblock.first;
          const int d1 = -(ispace->nholes() - source_space->nholes());
          const int d2 = ispace->nele2() - source_space->nele2();
          const int d3 = ispace->nparticles() - source_space->nparticles();

          // equivalent to checking that the three values are 0, +1, -1 in any order
          // i.e. creation and annihilation are in separate RASes
          if (d1 + d2 + d3 == 0 && d1*d1 + d2*d2 + d3*d3 == 2) {
            const int create = ( d1 == 1 ? 0 : (d2 == 1 ? 1 : 2) );
            const int annihilate = ( d1 == -1 ? 0 : (d2 == -1 ? 1 : 2) );
            const int other = 3 - (create + annihilate);

            shared_ptr<RASGraph>& left_create_graph = ispace->graph(create);
            shared_ptr<RASGraph>& right_create_graph = source_space->graph(create);

            // info corresponding to creation operator
            vector<tuple<size_t, size_t, int, bitset<nbit__>>> create_data;
            {
              size_t right_index = 0;
              for(auto& ibit : right_create_graph->model_space()) {
                for (int i = create_orbstart; i < create_orbend; ++i) {
                  if (!ibit[i]) {
                    bitset<nbit__> targetbit = ibit; targetbit.set(i);
                    const size_t left_index = left_create_graph.lexical(create_orbstart, create_orbend, targetbit);
                    create_data.emplace_back(left_index, right_index, i, targetbit);
                  }
                }
                ++right_index;
              }
              assert(!create_data.empty());
            }

            shared_ptr<RASGraph>& left_ann_graph = ispace->graph(annihilate);
            shared_ptr<RASGraph>& right_ann_graph = source_space->graph(annihilate);

            // annihilation operator
            vector<tuple<size_t, size_t, int, bitset<nbit__>>> annihilate_data;
            {
              size_t right_index = 0;
              for(auto& ibit : right_ann_graph->model_space()) {
                for (int i = ann_orbstart; i < ann_orbend; ++i) {
                  if (ibit[i]) {
                    bitset<nbit__> targetbit = ibit; targetbit.reset(i);
                    const size_t left_index = left_ann_graph.lexical(ann_orbstart, ann_orbend, targetbit);
                    annihilate_data.emplace_back(left_index, right_index, i, targetbit);
                  }
                }
                ++right_index;
              }
              assert(!annihilate_data.empty());
            }

            map<pair<int, int>, double>> coords;
            for (auto &iann : annihilate_data) {
              for (auto& icre : create_data) {
                const int col = get<1>(icre)*right_cstride + get<1>(iann)*right_astride;
                const int row = get<0>(icre)*left_cstride + get<0>(iann)*right_astride;
                const int k = get<2>(icre);
                const int l = get<2>(iann);
                const int sign = (create < annihilate ? det->sign(get<3>(icre) ^ (get<3>(iann) << create_norbs), k, l) :
                    det->sign(get<3>(iann) ^ (get<3>(icre) << ann_norbs), k, l) ) ;
                // making the transpose
                coords.emplace(make_pair(row, col), sign * mo2e_ij[norb*(k + norb*norb*l)]);
              }
            }

            const int nsize = ( create == 2 || annihilate == 2 ? /* full left space */ : /* truncated left space */ );
            const int msize = ( create == 2 || annihilate == 2 ? /* full left space */ : /* truncated left space */ );
            auto sparse = make_shared<SparseMatrix>(msize, nsize, coords);
            // TODO Everything above here could be done ahead of time and saved

            shared_ptr<Matrix> Cp = cpblock.second->trans();
            const int nmultiplies = ispace->graph(other)->size();
            for (int imult = 0; imult < nmultiplies; ++imult) {
              dcsrmm_("N", sparse->ndim(), Cp->mdim(), sparse->mdim(), 1.0, sparse->data(), sparse->cols(), sparse->rind(),
                  Cp->element_ptr(imult*ostride, imult*ostride), Cp->ndim(),
                  1.0, VI_trans->element_ptr(imult*ostride,0), VI_trans->ndim());
            }
          }
          // all of them are zero, i.e. 'diagonal'
          else if (d1*d1 + d2*d2 + d3*d3 == 0) {
            //build three dense matrices
            for(int iras = 0; iras < 3; ++iras) {
              shared_ptr<RASGraph>& subgraph = ispace->graph(iras);
              // TODO should consider holding onto these blocks somehow
              shared_ptr<Matrix> dense_block(subgraph->size(), subgraph->size());
              // TODO figure out orbstart/fence
              size_t source = 0;
              // TODO write modelspace() code
              for(auto& ibit : subgraph->modelspace()) {
                for (int k = orbstart; k < orbfence; ++k) {
                  if (!ibit[k]) continue;
                  bitset<nbit__> tmpbit = ibit; tmpbit.reset(k);
                  for (int l = orbstart; l < orbfence; ++l) {
                    if (tmpbit[l]) continue;
                    bitset<nbit__> targetbit = tmpbit; targetbit.set(l);
                    const size_t target = subgraph->lexical(orbstart, orbfence, targetbit);
                    const int sign = det->sign(ibit, k, l);
                    dense_block->element(target, source) += sign * mo2e[norb*(l + k*norb*norb)];
                  }
                }
                ++source;
              }

              // dense_block is made, now need to apply it
              shared_ptr<const Matrix> Cp = cpblock.second->get(iras);
              const size_t nblocks = cpblock.second->nblocks(iras);
              const size_t blocksize = dense_block->ndim();
              auto tmpV = make_shared<Matrix>(Cp->ndim(), blocksize);
              shared_ptr<Matrix> destination = VI[iras];
              for (size_t ib = 0; ib < nblocks; ++ib) {
                // verify that this is actually calling dgemm correctly
                dgemm_("N", "N", Cp->ndim(), dense_block->mdim(), dense_block->mdim(), 1.0, Cp->element_ptr(0, ib*blocksize), Cp->ndim(),
                    dense_block->data(), dense_block->ndim(), 0.0, tmpV->data(), tmpV->ndim());
                // add proper sublocks to destination
                for (auto& subblock : cpblock.second->dmap())
                  destination->add_strided_block(subblock.offset(), ib*blocksize, subblock.size(), blocksize, Cp->ndim(), tmpV->data());
              }
            }
          }
        }

        // TODO set nr1,2,3
        for (size_t r3 = 0, col = 0; r3 < nr3; ++r3) {
          for (size_t r1 = 0; r1 < nr1; ++r1) {
            for (size_t r2 = 0; r2 < nr2; ++r2, ++col) {
              double* targ = VI_out->element_ptr(0, col);
              blas::ax_plus_y_n(1.0, VI[0]->element_ptr(0, r1 + nr1*r2 + nr1*nr2*r3), phisize, targ);
              blas::ax_plus_y_n(1.0, VI[2]->element_ptr(0, r3 + nr3*r2 + nr3*nr2*r1), phisize, targ);
            }
          }
        }
        VI[0].reset(); VI[2].reset();

        // scatter
        double* vdata = VI_out->data();
        for (auto& iphiblock : det->phib_ij(ij) ) {
          for (auto& iphi : iphiblock) {
            shared_ptr<const StringSpace> betaspace = det->space<1>(det->stringb(iphi.target));
            if (det->allowed(ispace, betaspace)) {
              const double* sourcedata = vdata;

              shared_ptr<RASBlock<double>> sgblock = sigma->block(betaspace, ispace);
              double* targetdata = sgblock->data() + iphi.target - betaspace->offset();

              const size_t lb = sgblock->lenb();

              for (size_t i = 0; i < la; ++i) {
                targetdata[i*lb] += sourcedata[i*lb];
              }
            }

            ++vdata;
          }
        }
      }
    }
  }
}
#endif
