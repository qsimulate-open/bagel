//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest_prod_asd.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd_dmrg/gamma_forest_prod_asd.h>
#include <src/asd_dmrg/product_civec.h>
#include <src/math/btas_interface.h>

using namespace bagel;
using namespace std;

// collect information to get ready for compute, but does not run compute
GammaForestProdASD::GammaForestProdASD(map<BlockKey, vector<shared_ptr<const ProductRASCivec>>> block_states) : block_states_(block_states) {
  // set nele which will be used for the striding
  nele_ = 0;
  for (auto& i : block_states)
    nele_ = max(nele_, static_cast<size_t>(i.second.front()->nele()));

  // reorganize block_states into a map of Dvecs
  map<ProductState, shared_ptr<RASDvec>> vecmap;
  for (auto& state_block : block_states) {
    const int nstates = state_block.second.size();
    shared_ptr<const DMRG_Block> block = state_block.second.front()->left();
    for (auto& binfo : block->blocks()) {
      const int nb = binfo.nstates;
      for (int i = 0; i < nb; ++i) {
        vector<shared_ptr<RASCivec>> tmp;
        for (int j = 0; j < nstates; ++j)
          tmp.push_back(make_shared<RASCivec>(state_block.second.at(j)->sector(binfo.key())->civec(i)));
        auto dvec = make_shared<RASDvec>(tmp);
        vecmap.emplace(ProductState(binfo.key(), BlockKey(dvec->det()->nelea(), dvec->det()->neleb()), i), dvec);
      }
    }
  }

  possible_couplings_ = {
    {GammaSQ::CreateAlpha},
    {GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateAlpha,     GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateBeta,      GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateBeta,      GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateBeta,      GammaSQ::CreateAlpha,     GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateAlpha,     GammaSQ::CreateAlpha,     GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateBeta,      GammaSQ::CreateBeta,      GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateAlpha,     GammaSQ::CreateBeta,      GammaSQ::AnnihilateBeta},
  };

  forest_ = make_shared<GammaForest<RASDvec, 1>>();

  for (auto& ket : vecmap) {
    for (auto& coupling : possible_couplings_) {
      BlockKey new_key = apply_key(ket.first.key(), coupling);

      for (auto& bra : vecmap) {
        if (new_key==bra.first.key()) { // this is a POTENTIAL combination (although not guaranteed)
          // number of ways to partition the operators into either the block or the ci part
          const int npartition = 1 << coupling.size();

          for (int partition = 0; partition < npartition; ++partition) {
            // if bit is true --> ci part; if bit is false --> block part
            bitset<3> pattern(partition);
            list<GammaSQ> blockops, ciops;
            {
              auto citer = coupling.begin();
              for (int op = 0; op < coupling.size(); ++op, ++citer)
                (pattern[op] ? ciops : blockops).push_back(*citer);
            }
            if (apply_key(ket.first.block,blockops)==bra.first.block && apply_key(ket.first.ci,ciops)==bra.first.ci) {
              // this partitioning contributes, and now just need to figure out in what way to append it to the forest
              // block parts are in there already, so just need to worry about the ci part
              auto bra_ptr = &bra;
              auto ket_ptr = &ket;

              bool conj;
              list<GammaSQ> ordered_ciops;
              tie(conj, ignore, ordered_ciops) = try_permutations(ciops);

              if (!(ordered_ciops.size()==coupling.size() && bra_ptr->first.state!=ket_ptr->first.state)) {
                forest_->insert<0>(bra_ptr->second, state_tag(bra_ptr->first), ket_ptr->second, state_tag(ket_ptr->first), ordered_ciops);
              }
            }
          }
        }
      }
    }
  }
}

tuple</*conj*/bool, /*rev*/bool, list<GammaSQ>> GammaForestProdASD::try_permutations(const list<GammaSQ>& gammalist) const {
  if (gammalist.empty()) return make_tuple(false, false, gammalist);
  for (int conjrev = 0; conjrev < 4; ++conjrev) {
    const bool rev = bitset<2>(conjrev)[1];
    const bool conj = bitset<2>(conjrev)[0];

    list<GammaSQ> tmp = rev^conj ? list<GammaSQ>(gammalist.rbegin(), gammalist.rend()) : list<GammaSQ>(gammalist.begin(), gammalist.end());
    if (conj) for_each(tmp.begin(), tmp.end(), [] (GammaSQ& a) { a = conjugate(a); });

    if (count(possible_couplings_.begin(), possible_couplings_.end(), tmp)==1)
      return make_tuple(conj, rev, tmp);
  }
  assert(false);
  return make_tuple(false, false, list<GammaSQ>());
}

tuple<int, int, int> GammaForestProdASD::get_indices(const bitset<3> bit, const int size, const int ijk_local, const int lnorb, const bool block_is_reversed, const int rnorb, const bool ci_is_reversed) const {
  vector<int> strides(size);
  strides[0] = 1;
  for (int i = 1; i < size; ++i)
    strides[i] = strides[i-1] * (bit[i] ? rnorb : lnorb);

  vector<int> indices(size);
  int current = ijk_local;
  for (int i = size-1; i >= 0; --i) {
    const int index = current / strides[i];
    indices[i] = index;
    current -= index * strides[i];
  }
  assert(current==0);

  const int norb = lnorb + rnorb;

  const int ijk = accumulate(indices.rbegin(), indices.rend(), 1, [norb] (int ij, int index) { return index + ij*norb; });

  vector<int> block_indices, ci_indices;
  for (int i = 0; i < size; ++i) (bit[i] ? ci_indices : block_indices).push_back(indices[i]);
  if (block_is_reversed) reverse(block_indices.begin(), block_indices.end());
  if (ci_is_reversed) reverse(ci_indices.begin(), ci_indices.end());

  const int block_index = accumulate(block_indices.rbegin(), block_indices.rend(), 1, [lnorb] (int ij, int index) { return index + ij*lnorb; });
  const int ci_index = accumulate(ci_indices.rbegin(), ci_indices.rend(), 1, [rnorb] (int ij, int index) { return index + ij*rnorb; });

  return make_tuple(ijk, block_index, ci_index);
}

void GammaForestProdASD::compute() {
  // first hit compute on forest
  forest_->compute();

  const int lnorb = block_states_.begin()->second.front()->left()->norb();
  const int rnorb = block_states_.begin()->second.front()->space()->norb();
  const int norb = lnorb + rnorb;

  // loop over block states to allocate tree
  for (auto& ket_states : block_states_) {
    for (auto& coupling : possible_couplings_) {
      BlockKey bra_key = apply_key(ket_states.first, coupling);
      if (block_states_.find(bra_key)==block_states_.end()) continue;
      const int nijk = accumulate(coupling.begin(), coupling.end(), 1, [norb] (int x, GammaSQ a) { return x*norb; });
      const int nbrastates = block_states_.at(bra_key).size();
      const int nketstates = ket_states.second.size();
      auto gamma_matrix = make_shared<Matrix>(nbrastates*nketstates, nijk);

      shared_ptr<const DMRG_Block> dmrgblock = ket_states.second.front()->left();
      assert(dmrgblock == block_states_.at(bra_key).front()->left());

      for (auto& kbinfo : dmrgblock->blocks()) {
        // key for the RAS part of the wavefunction
        BlockKey krkey(ket_states.first.nelea - kbinfo.nelea, ket_states.first.neleb - kbinfo.neleb);

        // loop through all the ways to divide the operators between block and ras
        const int npart = 1 << coupling.size();
        for (int part = 0; part < npart; ++part) {
          bitset<3> bit(part);
          list<GammaSQ> blockops, ciops;
          auto citer = coupling.begin();
          // divide such that if bit is true --> ciops, bit false --> blockops
          for (int i = 0; i < coupling.size(); ++i, ++citer)
            (bit[i] ? ciops : blockops).push_back(*citer);

          // apply blockops to kbinfo
          BlockKey bbkey = apply_key(kbinfo.key(), blockops);
          if (dmrgblock->contains(bbkey)) {
            BlockInfo bbinfo = dmrgblock->blockinfo(bbkey);

            int Mbra = bbinfo.nstates;
            int Mket = kbinfo.nstates;

            // set up a bunch of key objects that will be used for the rest of this loop
            BlockKey block_bra = bbinfo.key();
            BlockKey block_ket = kbinfo.key();
            BlockKey ci_bra(bra_key.nelea - block_bra.nelea, bra_key.neleb - block_bra.neleb);
            BlockKey ci_ket(ket_states.first.nelea - block_ket.nelea, ket_states.first.neleb - block_ket.neleb);

            ProductState ps_bra(block_bra, ci_bra, 0);
            ProductState ps_ket(block_ket, ci_ket, 0);

            // figure out how the operations need to be reordered so they are in the set of computed couplings
            bool block_conj, block_rev, ci_conj, ci_rev;
            tie(block_conj, block_rev, blockops) = try_permutations(blockops);
            tie(ci_conj, ci_rev, ciops) = try_permutations(ciops);

            const int phase = (block_rev^ci_rev ? -1 : 1) * (part==2 || part==5 ? -1 : 1) * static_cast<int>(1 - (((blockops.size()*(block_bra.nelea+block_bra.neleb))%2) << 1));

            // swap where appropriate
            if (block_conj) swap(block_bra, block_ket);
            if (ci_conj) swap(ci_bra, ci_ket);
            if (ci_conj) swap(Mbra, Mket);
            if (ci_conj) swap(ps_bra, ps_ket);

            const int nijk_part = accumulate(blockops.begin(), blockops.end(), 1, [lnorb] (int x, GammaSQ a) { return x*lnorb; }) *
                                  accumulate(ciops.begin(), ciops.end(), 1, [rnorb] (int x, GammaSQ a) { return x*rnorb; });

            shared_ptr<const btas::Tensor3<double>> block_part = dmrgblock->coupling(blockops).at({block_bra, block_ket}).data;

            for (int ijk_part = 0; ijk_part < nijk_part; ++ijk_part) {
              int ijk, block_index, ci_index;
              tie(ijk, block_index, ci_index) = get_indices(bit, coupling.size(), ijk_part, lnorb, block_conj^block_rev, rnorb, ci_conj^ci_rev);

              // loop through all vectors in the sector
              for (int bra_k = 0; bra_k < Mbra; ++bra_k) {
                for (int ket_k = 0; ket_k < Mket; ++ket_k) {
                  if (blockops.empty() && ket_k!=bra_k) continue; // block states are orthonormal

                  ps_bra.state = bra_k;
                  ps_ket.state = ket_k;

                  shared_ptr<const Matrix> ras_part = forest_->get<0>(state_tag(ps_bra), state_tag(ps_ket), ciops);
                  const double block_gamma_ss = (*block_part)(bra_k, ket_k, block_index);

                  // fill in part of gamma
                  double* target = gamma_matrix->element_ptr(0, ijk);
                  if (ci_conj)
                    blas::transpose(ras_part->element_ptr(0, ci_index), nketstates, nbrastates, target, static_cast<double>(phase)*block_gamma_ss);
                  else
                    blas::ax_plus_y_n(static_cast<double>(phase)*block_gamma_ss, ras_part->element_ptr(0, ci_index), ras_part->ndim(), target);
                }
              }
            }
          }
        }
      }
      gammas_.emplace(make_tuple(coupling, bra_key, ket_states.first), gamma_matrix);
    }
  }
}
