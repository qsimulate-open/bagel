//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/form_sigma.cc
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

#include <map>
#include <src/util/math/sparsematrix.h>
#include <src/ci/ras/form_sigma.h>
#include <src/ci/ras/sparse_ij.h>

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

void FormSigmaRAS::sigma_ab(const RASCivecView cc, RASCivecView sigma, const double* mo2e) const {
  assert(*cc.det() == *sigma.det());
  shared_ptr<const RASDeterminants> det = cc.det();

  const int norb = det->norb();

  // pre-compute all sparse F matrices
  Sparse_IJ sparseij(det->stringspaceb(), det->stringspaceb());

  // figure out maximum block size
  const size_t max_ccblock_size = (*max_element(cc.blocks().begin(), cc.blocks().end(),
          [] (shared_ptr<const RASBlock<double>> a, shared_ptr<const RASBlock<double>> b) {
            return ( a ? a->size() : 0) < ( b ? b->size() : 0);
          }))->size();

  // allocate some scratch space. these upperbounds may be overkill
  unique_ptr<double[]> cprime(new double[2*max_ccblock_size]);
  unique_ptr<double[]> V(new double[2*max_ccblock_size]);

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
              const shared_ptr<const RASBlock<double>>& tblock = cc.block(target_bspace, target_aspace);
              const size_t o = tblock->offset() + (phi.target - target_aspace->offset()) * tblock->lenb();
              reduced_phi.emplace_back(phi.source, phi.sign, o);
            }
          }

          if (reduced_phi.empty()) continue;

          for (auto& source_block : cc.allowed_blocks<0>(source_aspace)) {
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
              fill_n(cprime.get(), slb * reduced_phi.size(), 0.0);
              int current = 0;
              for (auto& i : reduced_phi)
                blas::ax_plus_y_n(get<1>(i), source_block->data() + slb*get<0>(i), slb, cprime.get() + current++*slb);

              // compute V = F * C'
              dcsrmm_("N", tlb, reduced_phi.size(), slb, 1.0, sparseF->data(), sparseF->cols(), sparseF->rind(), cprime.get(), slb, 0.0, V.get(), tlb);

              // scatter to add V to sigma
              current = 0;
              for (auto& i : reduced_phi)
                blas::ax_plus_y_n(1.0, V.get() + tlb*current++, tlb, sigma.data() + get<2>(i));
            }
          }
        }
      }
    }
  }
}
