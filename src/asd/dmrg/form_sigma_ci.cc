//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/form_sigma.cc
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

#include <map>
#include <src/util/math/sparsematrix.h>
#include <src/asd/dmrg/form_sigma.h>
#include <src/ci/ras/form_sigma.h>
#include <src/ci/ras/apply_operator.h>

using namespace std;
using namespace bagel;

void FormSigmaProdRAS::resolve_H_aa(const RASBlockVectors& cc, RASBlockVectors& sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASDeterminants> det = cc.det();
  assert(*det == *sigma.det());

  const int norb = det->norb();
  const size_t la = det->lena();

  const int M = cc.mdim();
  assert(M == sigma.mdim());

#if HAVE_MPI_H
  StaticDist dist(M, mpi__->size());

  size_t mstart, mend;
  tie(mstart, mend) = dist.range(mpi__->rank());
#else
  const size_t mstart = 0;
  const size_t mend = M;
#endif

  Matrix F(la, min(static_cast<size_t>(batchsize_), la));

  for (auto& target_aspace : *det->stringspacea()) {
    // Do this multiplication batchwise
    const int nbatches = (target_aspace->size() - 1)/batchsize_ + 1;
    for (int batch = 0; batch < nbatches; ++batch) {
      const size_t batchstart = batch * batchsize_;
      const size_t batchlength = min(static_cast<size_t>(batchsize_), target_aspace->size() - batchstart);

      F.zero();

      for (size_t ia = 0; ia < batchlength; ++ia) {
        double * const fdata = F.element_ptr(0, ia);
        const size_t offset = batchstart + target_aspace->offset();
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
      for (auto& source_block : det->blockinfo()) {
        if (source_block->empty()) continue;
        if (!det->allowed(target_aspace, source_block->stringsb())) continue;
        const shared_ptr<const CIBlockInfo<RASString>>& target_block = det->blockinfo(source_block->stringsb(), target_aspace);

        assert(source_block->lenb() == target_block->lenb());
        for (int m = mstart; m < mend; ++m) {
          const double* sourcedata = cc.element_ptr(source_block->offset(), m);
          double* targetdata = sigma.element_ptr(target_block->offset() + batchstart*target_block->lenb(), m);
          dgemm_("N", "N", target_block->lenb(), batchlength, source_block->lena(), 1.0,
                           sourcedata, source_block->lenb(),
                           F.element_ptr(source_block->stringsa()->offset(), 0), F.ndim(), 1.0,
                           targetdata, target_block->lenb());
        }
      }
    }
  }
}


void FormSigmaProdRAS::resolve_H_bb(const RASBlockVectors& cc, RASBlockVectors& sigma, shared_ptr<const RASDeterminants> trans_det, const double* g, const double* mo2e) const {
  RASBlockVectors cc_trans = cc.transpose_civecs(trans_det);
  RASBlockVectors sigma_trans(cc_trans.det(), cc_trans.mdim());

  resolve_H_aa(cc_trans, sigma_trans, g, mo2e);

  sigma.ax_plus_y(1.0, sigma_trans.transpose_civecs(sigma.det()));
}



void FormSigmaProdRAS::resolve_H_ab(const RASBlockVectors& cc, RASBlockVectors& sigma, const Sparse_IJ& sparseij, const double* mo2e) const {
  assert(*cc.det() == *sigma.det());
  shared_ptr<const RASDeterminants> det = cc.det();

  const int norb = det->norb();
  const int M = cc.mdim();
  assert(M == sigma.mdim());

#if HAVE_MPI_H
  StaticDist dist(M, mpi__->size());

  size_t mstart, mend;
  tie(mstart, mend) = dist.range(mpi__->rank());
#else
  const size_t mstart = 0;
  const size_t mend = M;
#endif

  const size_t msize = mend - mstart;

  // figure out maximum block size
  const size_t max_ccblock_size = (*max_element(det->blockinfo().begin(), det->blockinfo().end(),
          [] (const shared_ptr<const CIBlockInfo<RASString>>& a, const shared_ptr<const CIBlockInfo<RASString>>& b) {
            return ( a ? a->size() : 0) < ( b ? b->size() : 0);
          }))->size();

  // allocate some scratch space. these upperbounds may be overkill
  unique_ptr<double[]> cprime(new double[2*max_ccblock_size*msize]);
  unique_ptr<double[]> V(new double[2*max_ccblock_size*msize]);

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      const double* mo2e_ij = mo2e + i + norb*norb*j;
      for (auto& target_bspace : *det->stringspaceb()) {
        const size_t tlb = target_bspace->size();
        // looping over source_aspace
        for (auto& phiblock : det->phia_ij(ij) ) {
          const shared_ptr<const RASString>& source_aspace = phiblock.source_space();

          // make a reduced list of only those excitations that will contribute to the sigma vector
          vector<tuple</*source*/size_t,/*sign*/int, /*offset_of_target*/size_t>> reduced_phi;
          for (auto& phi : phiblock) {
            auto target_aspace = det->space<0>(det->string_bits_a(phi.target));
            if (det->allowed(target_aspace, target_bspace)) {
              const shared_ptr<const CIBlockInfo<RASString>>& tblock = det->blockinfo(target_bspace, target_aspace);
              const size_t o = tblock->offset() + (phi.target - target_aspace->offset()) * tblock->lenb();
              reduced_phi.emplace_back(phi.source, phi.sign, o);
            }
          }

          if (reduced_phi.empty()) continue;

          for (auto& source_block : det->matching_blocks<0>(source_aspace)) {
            const shared_ptr<const RASString>& source_bspace = source_block->stringsb();
            const size_t slb = source_bspace->size();

            // F matrix in sparse format
            const shared_ptr<SparseMatrix>& sparseF = sparseij.sparse_matrix(target_bspace->tag(), source_bspace->tag());

            // if this assert fails, max_ccblock_size is not a good enough upper bound
            assert(max_ccblock_size >= max(slb, tlb) * reduced_phi.size());

            if (sparseF) {
              // fill in sparse matrix
              sparseF->zero();
              for (auto& iter : sparseij.sparse_data(target_bspace->tag(), source_bspace->tag()))
                *iter.ptr += static_cast<double>(iter.sign) * mo2e_ij[norb*(iter.i + norb*norb*iter.j)];

              // gather to fill in C'
              fill_n(cprime.get(), slb * reduced_phi.size() * msize, 0.0);
              int current = 0;
              for (int m = mstart; m < mend; ++m) {
                const double* sourcedata = cc.element_ptr(source_block->offset(), m);
                for (auto& i : reduced_phi)
                  blas::ax_plus_y_n(get<1>(i), sourcedata + slb*get<0>(i), slb, cprime.get() + current++*slb);
              }

              // compute V = F * C'
              dcsrmm_("N", tlb, reduced_phi.size()*msize, slb, 1.0, sparseF->data(), sparseF->cols(), sparseF->rind(), cprime.get(), slb, 0.0, V.get(), tlb);

              // scatter to add V to sigma
              current = 0;
              for (int m = mstart; m < mend; ++m)
                for (auto& i : reduced_phi)
                  blas::ax_plus_y_n(1.0, V.get() + tlb*current++, tlb, sigma.element_ptr(get<2>(i), m));
            }
          }
        }
      }
    }
  }
}


void FormSigmaProdRAS::resolve_S_aaa(const RASBlockVectors& cc, RASBlockVectors& sigma, const double* Jp, const PhiIJKLists& phi_ijk) const {
  shared_ptr<const RASDeterminants> sdet = cc.det();
  shared_ptr<const RASDeterminants> tdet = sigma.det();

  const int norb = sdet->norb();
  assert(norb == tdet->norb());

  const size_t batchsize = batchsize_;
  const size_t sla = sdet->lena();

  const int M = cc.mdim();
  assert(M == sigma.mdim());

  Matrix F(sla, min(batchsize, tdet->lena()), true);
  for (auto& targetspace : *tdet->stringspacea()) {
    const int nbatches = (targetspace->size()-1)/batchsize + 1;
    for (int batch = 0; batch < nbatches; ++batch) {
      const size_t batchstart = batch * batchsize;
      const size_t batchlength = min(batchsize, targetspace->size() - batchstart);

      F.zero();
      for (size_t ia = 0; ia < batchlength; ++ia) {
        double* const fdata = F.element_ptr(0, ia);
        for (auto& iter : phi_ijk.data(ia + batchstart + targetspace->offset()))
          fdata[iter.source] += static_cast<double>(iter.sign) *
                      (Jp[iter.j+norb*iter.i+norb*norb*iter.k] - Jp[iter.k+norb*iter.i+norb*norb*iter.j]);
      }

      for (auto& ccblock : sdet->blockinfo()) {
        if (ccblock->empty()) continue;
        if (!tdet->allowed(targetspace, ccblock->stringsb())) continue;
        const shared_ptr<const CIBlockInfo<RASString>>& target_block = tdet->blockinfo(ccblock->stringsb(), targetspace);

        assert(ccblock->lenb() == target_block->lenb());
        const double* fdata = F.element_ptr(ccblock->stringsa()->offset(), 0);
        for (int m = 0; m < M; ++m) {
          const double* source_data = cc.element_ptr(ccblock->offset(), m);
          double* target_data = sigma.element_ptr(target_block->offset() + batchstart * target_block->lenb(), m);
          dgemm_("N", "N", target_block->lenb(), batchlength, ccblock->lena(), 1.0, source_data, ccblock->lenb(),
                    fdata, F.ndim(), 1.0, target_data, target_block->lenb());
        }
      }
    }
  }
}


// computes \sum_{ijk} (k)_alpha i^+_beta j_beta (pk|ij)
// (k) can be creation or annihilation and is determined based on the electron count of sigma
void FormSigmaProdRAS::resolve_S_abb(const RASBlockVectors& cc, RASBlockVectors& sigma, const double* Jp, const PhiKLists& phik, const Sparse_IJ& sparseij) const {
  shared_ptr<const RASDeterminants> sdet = cc.det();
  shared_ptr<const RASDeterminants> tdet = sigma.det();

  assert(abs(sdet->nelea()-tdet->nelea())==1);
  assert(sdet->neleb()==tdet->neleb());

  const int M = cc.mdim();
  assert(M == sigma.mdim());

  // first, figure out maximum size of all the blocks
  const size_t max_ccblock_size = (*max_element(sdet->blockinfo().begin(), sdet->blockinfo().end(),
          [] (const shared_ptr<const CIBlockInfo<RASString>>& a, const shared_ptr<const CIBlockInfo<RASString>>& b) {
            return ( a ? a->size() : 0) < ( b ? b->size() : 0);
          }))->size() * M;
  const size_t max_sgblock_size = (*max_element(tdet->blockinfo().begin(), tdet->blockinfo().end(),
          [] (const shared_ptr<const CIBlockInfo<RASString>>& a, const shared_ptr<const CIBlockInfo<RASString>>& b) {
            return ( a ? a->size() : 0) < ( b ? b->size() : 0);
          }))->size() * M;

  // allocate chunks of storage equal to the maximum possible size that may be needed. probably overkill, but also probably fine
  unique_ptr<double[]> cprime(new double[max_ccblock_size]);
  unique_ptr<double[]> V(new double[max_sgblock_size]);

  const int norb = sdet->norb();
  assert(norb == tdet->norb());

  // k^?_alpha i^+_beta j_beta portion. the harder part
  for (int k = 0; k < norb; ++k) {
    for (auto& target_bspace : *tdet->stringspaceb()) {
      const size_t tlb = target_bspace->size();
      for (auto& source_aspace : *sdet->stringspacea()) {
        const vector<PhiKLists::PhiK>& full_phi = phik.data(k).at(source_aspace->tag());
        vector<PhiKLists::PhiK> reduced_RI;
        reduced_RI.reserve(count_if(full_phi.begin(), full_phi.end(), [&tdet, &target_bspace] (const PhiKLists::PhiK& i) {
          const RASString* target_aspace = i.target_space;
          return tdet->allowed(target_aspace->nholes(), target_bspace->nholes(), target_aspace->nparticles(), target_bspace->nparticles());
        }));

        // after this step, the "target" field of the PhiK objects gives the starting position within the civector
        for (auto& i : full_phi) {
          const RASString* target_aspace = i.target_space;
          if (tdet->allowed(target_aspace->nholes(), target_bspace->nholes(), target_aspace->nparticles(), target_bspace->nparticles())) {
            const shared_ptr<const CIBlockInfo<RASString>>& bi = tdet->blockinfo(target_aspace->nholes(), target_bspace->nholes(),
                                                                                 target_aspace->nparticles(), target_bspace->nparticles());
            reduced_RI.emplace_back(i.source, bi->offset() + i.target*bi->lenb(), i.sign, target_aspace);
          }
        }

        if (reduced_RI.empty()) continue;

        for (auto& source_block : sdet->matching_blocks<0>(source_aspace)) {
          auto& source_bspace = source_block->stringsb();
          const size_t slb = source_bspace->size();

          // Now build an F matrix in sparse format
          const shared_ptr<SparseMatrix>& sparseF = sparseij.sparse_matrix(target_bspace->tag(), source_bspace->tag());

          // if this assert fails, max_ccblock_size is not a good enough upperbound
          assert(max_ccblock_size >= slb * reduced_RI.size() * M);

          // if this assert fails, max_sgblock_size is not a good enough upperbound
          assert(max_sgblock_size >= tlb * reduced_RI.size() * M);

          if (sparseF) {
            sparseF->zero();
            for (auto& iter : sparseij.sparse_data(target_bspace->tag(), source_bspace->tag()))
              *iter.ptr += static_cast<double>(iter.sign) * Jp[iter.j + norb*iter.i + norb*norb*k];

            fill_n(cprime.get(), slb * reduced_RI.size() * M, 0.0);

            int current = 0;
            for (int m = 0; m < M; ++m) {
              const double* sourcedata = cc.element_ptr(source_block->offset(), m);

              for (auto& i : reduced_RI)
                blas::ax_plus_y_n(i.sign, sourcedata + slb*i.source, slb, cprime.get() + current++*slb);
            }

            dcsrmm_("N", tlb, reduced_RI.size() * M, slb, 1.0, sparseF->data(), sparseF->cols(), sparseF->rind(), cprime.get(), slb, 0.0, V.get(), tlb);

              // scatter
              current = 0;
            for (int m = 0; m < M; ++m) {
              for (auto& i : reduced_RI) {
                double* targetdata = sigma.element_ptr(i.target, m);
                blas::ax_plus_y_n(1.0, V.get() + tlb*current++, tlb, targetdata);
              }
            }
          }
        }
      }
    }
  }
}
