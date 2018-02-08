//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_rdmderiv.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/ci/fci/fci.h>
#include <src/wfn/rdm.h>
#include <src/util/math/algo.h>
#include <src/util/parallel/rmawindow.h>
#include <src/util/taskqueue.h>

using namespace std;
using namespace bagel;


shared_ptr<Dvec> FCI::rdm1deriv(const int target) const {

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(target);

  // 1RDM ci derivative
  // <I|E_ij|0>

  auto dbra = make_shared<Dvec>(cbra->det(), norb_ * norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  return dbra;
}


shared_ptr<Dvec> FCI::rdm2deriv(const int target) const {
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(target);

  // make  <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_ * norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // second make <J|E_kl|I><I|E_ij|0> - delta_li <J|E_kj|0>
  auto ebra = make_shared<Dvec>(cbra->det(), norb_ * norb_ * norb_ * norb_);
  auto tmp = make_shared<Dvec>(cbra->det(), norb_ * norb_);
  int ijkl = 0;
  int ij = 0;
  for (auto iter = dbra->dvec().begin(); iter != dbra->dvec().end(); ++iter, ++ij) {
    const int j = ij / norb_;
    const int i = ij - j*norb_;
    tmp->zero();
    sigma_2a1(*iter, tmp);
    sigma_2a2(*iter, tmp);
    int kl = 0;
    for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
      *ebra->data(ijkl) = **t;
      const int l = kl / norb_;
      const int k = kl - l * norb_;
      if (l == i) *ebra->data(ijkl) -= *dbra->data(k + norb_ * j);
    }
  }
  return ebra;
}


namespace bagel {
  class RDM2derivTask {
    protected:
      const int ij_;
      const int kl_;
      shared_ptr<Matrix> emat_;
      shared_ptr<const Matrix> dmat_;
      shared_ptr<const Dvec> cc_;
      const int norb_;
      const int size_;
      const int offset_;
    public:
      RDM2derivTask(const int ij, const int kl, shared_ptr<Matrix> e, shared_ptr<const Matrix> d, shared_ptr<const Dvec> cc, const int norb, const int size, const int offset)
        : ij_(ij), kl_(kl), emat_(e), dmat_(d), cc_(cc), norb_(norb), size_(size), offset_(offset) { }

      void compute() {
        const int lena = cc_->det()->lena();
        const int lenb = cc_->det()->lenb();
        const int norb2 = norb_ * norb_;
        const int j = ij_ / norb_;
        const int i = ij_ - j * norb_;
        const int l = kl_ / norb_;
        const int k = kl_ - l * norb_;
        const int klij = kl_ + ij_ * norb2;

        for (auto& iter : cc_->det()->phia(k, l)) {
          const size_t iaJ = iter.source;
          const size_t iaI = iter.target;
          const double sign = static_cast<double>(iter.sign);
          for (size_t ib = 0; ib != lenb; ++ib) {
            const size_t iJ = ib + iaJ * lenb;
            const size_t iI = ib + iaI * lenb;
            if ((iJ - offset_) < size_ && iJ >= offset_)
              emat_->element(iJ-offset_, klij) += sign * dmat_->element(iI, ij_);
          }
        }

        for (size_t ia = 0; ia != lena; ++ia) {
          for (auto& iter : cc_->det()->phib(k, l)) {
            const size_t ibJ = iter.source;
            const size_t ibI = iter.target;
            const double sign = static_cast<double>(iter.sign);
            const size_t iJ = ibJ + ia * lenb;
            const size_t iI = ibI + ia * lenb;
            if ((iJ - offset_) < size_ && iJ >= offset_)
              emat_->element(iJ-offset_, klij) += sign * dmat_->element(iI, ij_);
          }
        }

        if (i == l) {
          const int kj = k + j * norb_;
          blas::ax_plus_y_n(-1.0, dmat_->element_ptr(offset_, kj), size_, emat_->element_ptr(0, klij));
        }
      }
  };
}


shared_ptr<Matrix> FCI::rdm2deriv_offset(const int target, const size_t offset, const size_t size, shared_ptr<const Matrix> dmat, const bool parallel) const {
  const size_t norb2 = norb_ * norb_;
  auto emat = make_shared<Matrix>(size, norb2*norb2, /*local=*/!parallel);

  TaskQueue<RDM2derivTask> task(norb2*(norb2+1)/2);
  for (int ij = 0; ij != norb2; ++ij) {
    for (int kl = ij; kl != norb2; ++kl) {
      if (((kl-ij) % mpi__->size() == mpi__->rank()) || !parallel)
        task.emplace_back(ij, kl, emat, dmat, cc_, norb_, size, offset);
    }
  }
  task.compute();

  if (parallel)
    emat->allreduce();

  for (size_t ij = 0; ij != norb2; ++ij)
    for (size_t kl = 0; kl != ij; ++kl) {
      const size_t klij = kl + ij * norb2;
      const size_t ijkl = ij + kl * norb2;
      for (size_t iJ = 0; iJ != size; ++iJ)
        emat->element(iJ, klij) = emat->element(iJ, ijkl);
    }

  return emat;
}


namespace bagel {
  class RDM3derivTask {
    protected:
      const int ij_;
      const int kl_;
      shared_ptr<Matrix> fock_fbra_;
      shared_ptr<const Matrix> fock_ebra_mat_;
      shared_ptr<const Matrix> ebra_;
      shared_ptr<const Matrix> fock_;
      shared_ptr<const Dvec> cc_;
      const int norb_;
      const int size_;
      const int offset_;
    public:
      RDM3derivTask(const int ij, const int kl, shared_ptr<Matrix> f, shared_ptr<const Matrix> fe, shared_ptr<const Matrix> e, shared_ptr<const Matrix> fock, shared_ptr<const Dvec> cc, const int norb, const int size, const int offset)
        : ij_(ij), kl_(kl), fock_fbra_(f), fock_ebra_mat_(fe), ebra_(e), fock_(fock), cc_(cc), norb_(norb), size_(size), offset_(offset) { }

      void compute() {
        const size_t norb2 = norb_ * norb_;
        const size_t norb3 = norb2 * norb_;
        const int lena = cc_->det()->lena();
        const int lenb = cc_->det()->lenb();
        const int j = ij_ / norb_;
        const int i = ij_ - j * norb_;
        const int l = kl_ / norb_;
        const int k = kl_ - l * norb_;
        const int klij = kl_ + ij_ * norb2;

        for (auto& iter : cc_->det()->phia(k, l)) {
          const size_t iaJ = iter.source;
          const size_t iaI = iter.target;
          const double sign = static_cast<double>(iter.sign);

          for (size_t ib = 0; ib != lenb; ++ib) {
            const size_t iI = ib + iaI * lenb;
            const size_t iJ = ib + iaJ * lenb;
            if ((iJ - offset_) < size_ && iJ >= offset_)
              fock_fbra_->element(iJ-offset_, klij) += sign * fock_ebra_mat_->element(iI, ij_);
          }
        }

        for (size_t ia = 0; ia != lena; ++ia) {
          for (auto& iter : cc_->det()->phib(k, l)) {
            const size_t ibJ = iter.source;
            const size_t ibI = iter.target;
            const double sign = static_cast<double>(iter.sign);
            const size_t iI = ibI + ia * lenb;
            const size_t iJ = ibJ + ia * lenb;
            if ((iJ - offset_) < size_ && iJ >= offset_)
              fock_fbra_->element(iJ-offset_, klij) += sign * fock_ebra_mat_->element(iI, ij_);
          }
        }

        if (i == l) {
          const int kj = k + j * norb_;
          blas::ax_plus_y_n(-1.0, fock_ebra_mat_->element_ptr(offset_, kj), size_, fock_fbra_->element_ptr(0, klij));
        }

        for (int n = 0; n != norb_; ++n) {
          const size_t ijkn = j * norb3 + i * norb2 + n * norb_ + k;
          for (size_t iJ = 0; iJ != size_; ++iJ)
            fock_fbra_->element(iJ, klij) -= ebra_->element(iJ, ijkn) * fock_->element(l, n);
        }
      }
  };
}


shared_ptr<Matrix> FCI::rdm2fderiv(const int target, shared_ptr<const Matrix> fock, shared_ptr<const Matrix> dmat) const {
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(target);

  const size_t norb2 = norb_ * norb_;
  const size_t norb4 = norb2 * norb2;

  // Fock-weighted 2RDM derivative construction is multipassed and parallelized:
  //  (1) When ndet > 10000, (ndet < 10000 -> so fast, almost no gain)
  //  and (2) When we have processes more than one
  //  OR  (3) When the number of words in <I|E_ij,kl|0> is larger than (10,10) case (635,040,000)
  const size_t ndet = cbra->det()->size();
  const size_t ijmax = 635040001;
  const size_t ijnum = ndet * norb4;
  const size_t npass = ((mpi__->size() * 2 > ((ijnum - 1)/ijmax + 1)) && (mpi__->size() != 1) && ndet > 10000) ? mpi__->size() * 2 : (ijnum - 1) / ijmax + 1;
  const size_t nsize = (ndet - 1) / npass + 1;

  // form [J|k+l|0] = <J|m+k+ln|0> f_mn (multipassing with <J| )
  auto fock_emat = make_shared<Matrix>(ndet, norb2);
  for (size_t ipass = 0; ipass != npass; ++ipass) {
    if (ipass % mpi__->size() != mpi__->rank()) continue;

    const size_t ioffset = ipass * nsize;
    const size_t isize = (ipass != (npass - 1)) ? nsize : ndet - ioffset;
    shared_ptr<Matrix> ebra;
    ebra = rdm2deriv_offset(target, ioffset, isize, dmat, /*parallel=*/false);
    for (size_t kl = 0; kl != norb2; ++kl) {
      for (size_t mn = 0; mn != norb2; ++mn) {
        const size_t n = mn / norb_;
        const size_t m = mn - n * norb_;
        const size_t klmn = kl * norb2 + mn;
        for (size_t iI = 0; iI != isize; ++iI)
          fock_emat->element(iI + ioffset, kl) += fock->element(m, n) * ebra->element(iI, klmn);
      }
    }
  }
  fock_emat->allreduce();

  return fock_emat;
}


tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>
  FCI::rdm3deriv(const int target, shared_ptr<const Matrix> fock, const size_t offset, const size_t size, shared_ptr<const Matrix> dmat, shared_ptr<const Matrix> fock_emat) const {
#ifndef HAVE_MPI_H
  throw logic_error("FCI::rdm3deriv should not be called without MPI");
#endif
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(target);

  const size_t norb2 = norb_ * norb_;
  const size_t norb4 = norb2 * norb2;

  // then form [L|k+i+jl|0] <- <L|i+j|K>[K|k+l|0] + ...
  auto fock_fbra = make_shared<Matrix>(size, norb4, /*local=*/false);
  // Set ebra within the pass. This is 2RDM derivative
  auto ebra = rdm2deriv_offset(target, offset, size, dmat);

  // RDM3deriv contruction is task-base threaded
  TaskQueue<RDM3derivTask> task(norb2 * (norb2 + 1) / 2);
  for (int ij = 0; ij != norb2; ++ij) {
    for (int kl = ij; kl != norb2; ++kl) {
      if ((kl-ij) % mpi__->size() == mpi__->rank())
        task.emplace_back(ij, kl, fock_fbra, fock_emat, ebra, fock, cc_, norb_, size, offset);
    }
  }
  task.compute();
  fock_fbra->allreduce();

  for (size_t ij = 0; ij != norb2; ++ij)
    for (size_t kl = 0; kl != ij; ++kl) {
      const size_t klij = kl + ij * norb2;
      const size_t ijkl = ij + kl * norb2;
      for (size_t iJ = 0; iJ != size; ++iJ)
        fock_fbra->element(iJ, klij) = fock_fbra->element(iJ, ijkl);
    }

  return make_tuple(ebra, fock_fbra);
}
