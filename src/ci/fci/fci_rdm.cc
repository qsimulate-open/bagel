//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_rdm.cc
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
#include <src/util/prim_op.h>
#include <src/util/math/algo.h>
#include <src/wfn/rdm.h>

using namespace std;
using namespace bagel;

FCI_bare::FCI_bare(shared_ptr<const CIWfn> ci) {
  print_thresh_ = 1.0e-8;
  nelea_ = ci->det()->nelea();
  neleb_ = ci->det()->neleb();
  ncore_ = ci->ncore();
  norb_  = ci->nact();
  nstate_ = ci->nstates();
  energy_ = ci->energies();
  cc_ = ci->civectors() ? ci->civectors()->copy() : nullptr;
  det_ = ci->det();
  rdm1_ = make_shared<VecRDM<1>>();
  rdm2_ = make_shared<VecRDM<2>>();
}


void FCI::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_->at(0).
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norb_);
    rdm2_av_ = make_shared<RDM<2>>(norb_);
  } else if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  // we need expanded lists
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, /*compressed=*/false, /*mute=*/true);
  cc_->set_det(detex);

  for (int i = 0; i != nstate_; ++i)
    compute_rdm12(i, i);

  // calculate state averaged RDMs
  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weight_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weight_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0,0);
    rdm2_av_ = rdm2_->at(0,0);
  }

  cc_->set_det(det_);
}


void FCI::compute_rdm12(const int ist, const int jst) {
  if (det_->compress()) {
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }

  shared_ptr<Civec> ccbra = cc_->data(ist);
  shared_ptr<Civec> ccket = cc_->data(jst);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_civec(ccbra, ccket);

  // setting to private members.
  rdm1_->emplace(ist, jst, rdm1);
  rdm2_->emplace(ist, jst, rdm2);

  cc_->set_det(det_);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  FCI::compute_rdm12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {

  // since we consider here number conserving operators...
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb_*norb_);
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  return compute_rdm12_last_step(dbra, dket, cbra);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
FCI::compute_rdm12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {

  const int nri = cibra->asize()*cibra->lenb();

  // 1RDM c^dagger <I|\hat{E}|0>
  // 2RDM \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dketv = btas::group(*dket,0,2);
    auto rdm1t = btas::group(*rdm1,0,2);
    btas::contract(1.0, dketv, {0,1}, *cibra_data, {0}, 0.0, rdm1t, {1});

    auto dbrav = btas::group(*dbra,0,2);
    auto rdm2t = group(group(*rdm2, 2,4), 0,2);
    btas::contract(1.0, dbrav, {1,0}, dketv, {1,2}, 0.0, rdm2t, {0,2});
  }

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, rdm2->element_ptr(0,0,k,i));
    }
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb_; ++i)
    for (int k = 0; k != norb_; ++k)
      for (int j = 0; j != norb_; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
FCI::compute_rdm12_av_from_dvec(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Determinants> o) const {

  if (o != nullptr) {
    dbra->set_det(o);
    dket->set_det(o);
  }

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);

  assert(dbra->ij() == dket->ij() && dbra->det() == dket->det());

  for (int i = 0; i != dbra->ij(); ++i) {
    shared_ptr<RDM<1>> r1;
    shared_ptr<RDM<2>> r2;
    tie(r1, r2) = compute_rdm12_from_civec(dbra->data(i), dket->data(i));
    rdm1->ax_plus_y(weight_[i], r1);
    rdm2->ax_plus_y(weight_[i], r2);
  }

  if (o != nullptr) {
    dbra->set_det(det_);
    dket->set_det(det_);
  }

  return tie(rdm1, rdm2);
}

// computes 3 and 4RDM
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> FCI::rdm34(const int ist, const int jst) const {
  auto rdm3 = make_shared<RDM<3>>(norb_);
  auto rdm4 = make_shared<RDM<4>>(norb_);

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);

  shared_ptr<Civec> cbra = cc_->data(ist);
  shared_ptr<Civec> cket = cc_->data(jst);

  // first make <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket = dbra;
  if (cbra != cket) {
    dket = dbra->clone();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  }

  // second make <J|E_kl|I><I|E_ij|0> - delta_il <J|E_kj|0>
  auto make_evec_half = [this](shared_ptr<Dvec> d, shared_ptr<Matrix> e, const int dsize, const int offset) {
    const int norb2 = norb_ * norb_;
    const int lena = cc_->det()->lena();
    const int lenb = cc_->det()->lenb();
    int no = 0;

    for (int ij = 0; ij != norb2; ++ij) {
      const int j = ij/norb_;
      const int i = ij-j*norb_;

      for (int kl = ij; kl != norb2; ++kl) {
        const int l = kl/norb_;
        const int k = kl-l*norb_;

        for (auto& iter : cc_->det()->phia(k,l)) {
          size_t iaJ = iter.source;
          size_t iaI = iter.target;
          double sign = static_cast<double>(iter.sign);
          for (size_t ib = 0; ib != lenb; ++ib) {
            size_t iI = ib + iaI*lenb;
            size_t iJ = ib + iaJ*lenb;
            if ((iJ - offset) < dsize && iJ >= offset)
              e->element(iJ-offset, no) += sign * d->data(ij)->data(iI);
          }
        }

        for (size_t ia = 0; ia != lena; ++ia) {
          for (auto& iter : cc_->det()->phib(k,l)) {
            size_t ibJ = iter.source;
            size_t ibI = iter.target;
            double sign = static_cast<double>(iter.sign);
            size_t iI = ibI + ia*lenb;
            size_t iJ = ibJ + ia*lenb;
            if ((iJ - offset) < dsize && iJ >= offset)
              e->element(iJ-offset, no) += sign * d->data(ij)->data(iI);
          }
        }

        if (i == l) {
          const int kj = k+j*norb_;
          for (size_t iJ = offset; iJ != offset+dsize; ++iJ) {
            e->element(iJ-offset, no) -= d->data(kj)->data(iJ);
          }
        }
        ++no;
      }
    }
  };

  // RDM3, RDM4 construction is multipassed and parallelized:
  //  (1) When ndet > 10000, (ndet < 10000 -> too small, almost no gain)
  //  and (2) When we have processes more than one
  //  OR  (3) When the number of words in <I|E_ij,kl|0> is larger than (10,10) case (635,040,000)
  const size_t ndet = cbra->det()->size();
  const size_t norb2 = norb_ * norb_;
  const size_t ijmax = 635040001 * 2;
  const size_t ijnum = ndet * norb2 * norb2;
  const size_t npass = ((mpi__->size() * 2 > ((ijnum-1)/ijmax + 1)) && (mpi__->size() != 1) && ndet > 10000) ? mpi__->size() * 2 : (ijnum-1) / ijmax + 1;
  const size_t nsize = (ndet-1) / npass + 1;
  Timer timer;
  if (npass > 1) {
    cout << "    * Third and fourth order RDM (" << setw(2) << ist + 1 << "," << setw(2) << jst + 1 << ") evaluation" << endl;
    cout << "      will be done with " << npass << " passes" << endl;
  }

  rdm3->zero();
  rdm4->zero();

  // multipassing through {I}
  for (size_t ipass = 0; ipass != npass; ++ipass) {
    if (ipass % mpi__->size() != mpi__->rank()) continue;

    const size_t ioffset = ipass * nsize;
    const size_t isize = (ipass != (npass - 1)) ? nsize : ndet - ioffset;
    const size_t halfsize = norb2 * (norb2 + 1) / 2;
    auto eket_half = make_shared<Matrix>(isize, halfsize, /*local=*/true);
    make_evec_half(dket, eket_half, isize, ioffset);

    auto dbram = make_shared<Matrix>(isize, norb2, /*local=*/true);
    for (size_t ij = 0; ij != norb2; ++ij)
      copy_n(&(dbra->data(ij)->data(ioffset)), isize, dbram->element_ptr(0, ij));

    // put in third-order RDM: <0|E_mn|I><I|E_ij < kl|0>
    auto tmp3 = make_shared<Matrix>(*dbram % *eket_half);
    auto tmp3_full = make_shared<Matrix>(norb2, norb2 * norb2, /*local=*/true);

    for (size_t mn = 0; mn != norb2; ++mn) {
      int no = 0;
      for (size_t ij = 0; ij != norb2; ++ij) {
        for (size_t kl = ij; kl != norb2; ++kl) {
          int ijkl = ij + kl*norb2;
          int klij = kl + ij*norb2;
          tmp3_full->element(mn, ijkl) = tmp3->element(mn, no);
          tmp3_full->element(mn, klij) = tmp3->element(mn, no);
          ++no;
        }
      }
    }
    sort_indices<1,0,2,1,1,1,1>(tmp3_full->data(), rdm3->data(), norb_, norb_, norb2*norb2);

    // put in fourth-order RDM: <0|E_ij > kl|I><I|E_mn < op|0>
    shared_ptr<Matrix> ebra_half = eket_half;
    if (cbra != cket) {
      ebra_half = eket_half->clone();
      make_evec_half(dbra, ebra_half, isize, ioffset);
    }
    auto tmp4 = make_shared<Matrix>(*ebra_half % *eket_half);
    auto tmp4_full = make_shared<Matrix>(norb2 * norb2, norb2 * norb2, /*local=*/true);

    {
      int nklij = 0;
      for (size_t kl = 0; kl != norb2; ++kl) {
        for (size_t ij = kl; ij != norb2; ++ij) {
          int nmnop = 0;
          int klij = kl+ij*norb2;
          int ijkl = ij+kl*norb2;
          for (size_t mn = 0; mn != norb2; ++mn) {
            for (size_t op = mn; op != norb2; ++op) {
              int mnop = mn + op*norb2;
              int opmn = op + mn*norb2;
              tmp4_full->element(ijkl, mnop) = tmp4->element(nklij, nmnop);
              tmp4_full->element(klij, mnop) = tmp4->element(nklij, nmnop);
              tmp4_full->element(ijkl, opmn) = tmp4->element(nklij, nmnop);
              tmp4_full->element(klij, opmn) = tmp4->element(nklij, nmnop);
              ++nmnop;
            }
          }
          ++nklij;
        }
      }
    }
    sort_indices<1,0,3,2,4,1,1,1,1>(tmp4_full->data(), rdm4->data(), norb_, norb_, norb_, norb_, norb2*norb2);
  }

  rdm3->allreduce();
  rdm4->allreduce();

  if (npass > 1) {
    timer.tick_print("RDM evaluation (multipassing)");
  }

  // The remaining terms can be evaluated without multipassing
  {
    // then perform Eq. 49 of JCP 89 5803 (Werner's MRCI paper)
    // we assume that rdm2_[ist] is set
    for (int i0 = 0; i0 != norb_; ++i0)
      for (int i1 = 0; i1 != norb_; ++i1)
        for (int i2 = 0; i2 != norb_; ++i2)
          for (int i3 = 0; i3 != norb_; ++i3) {
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i2, i1, i0), norb_, rdm3->element_ptr(0, i3, i3, i2, i1, i0));
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i0, i3, i2), norb_, rdm3->element_ptr(0, i1, i3, i2, i1, i0));
          }
  }

  {
    for (int l = 0; l != norb_; ++l)
      for (int k = 0; k != norb_; ++k)
        for (int j = 0; j != norb_; ++j)
          for (int b = 0; b != norb_; ++b) {
            blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,0,0,k,b,l), norb_*norb_*norb_, rdm4->element_ptr(0,0,0,j,j,k,b,l));
            blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,0,0,l,b,k), norb_*norb_*norb_, rdm4->element_ptr(0,0,0,j,b,k,j,l));
            for (int i = 0; i != norb_; ++i) {
              blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0,k,b,l), norb_, rdm4->element_ptr(0,i,b,j,i,k,j,l));
              blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0,l,b,k), norb_, rdm4->element_ptr(0,i,b,j,j,k,i,l));
              for (int d = 0; d != norb_; ++d) {
                blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,k,b,j,d,l), norb_, rdm4->element_ptr(0,i,b,j,i,k,d,l));
                blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,l,b,j,d,k), norb_, rdm4->element_ptr(0,i,b,j,d,k,i,l));
              }
            }
          }
  }

  cc_->set_det(det_);

  return make_tuple(rdm3, rdm4);
}

// computes 3 and Fock-weighted 4RDM
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<3>>> FCI::rdm34f(const int ist, const int jst, shared_ptr<const Matrix> fock) const {
  auto rdm3 = make_shared<RDM<3>>(norb_);
  auto rdm4f = make_shared<RDM<3>>(norb_);

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);

  shared_ptr<Civec> cbra = cc_->data(ist);
  shared_ptr<Civec> cket = cc_->data(jst);

  // first make <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket = dbra;
  if (cbra != cket) {
    dket = dbra->clone();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  }

  // second make <J|E_kl|I><I|E_ij|0> - delta_il <J|E_kj|0>
  auto make_evec_half = [this](shared_ptr<Dvec> d, shared_ptr<Matrix> e, const int dsize, const int offset) {
    const int norb2 = norb_ * norb_;
    const int lena = cc_->det()->lena();
    const int lenb = cc_->det()->lenb();
    int no = 0;

    for (int ij = 0; ij != norb2; ++ij) {
      const int j = ij/norb_;
      const int i = ij-j*norb_;

      for (int kl = ij; kl != norb2; ++kl) {
        const int l = kl/norb_;
        const int k = kl-l*norb_;

        for (auto& iter : cc_->det()->phia(k,l)) {
          size_t iaJ = iter.source;
          size_t iaI = iter.target;
          double sign = static_cast<double>(iter.sign);
          for (size_t ib = 0; ib != lenb; ++ib) {
            size_t iI = ib + iaI*lenb;
            size_t iJ = ib + iaJ*lenb;
            if ((iJ - offset) < dsize && iJ >= offset)
              e->element(iJ-offset, no) += sign * d->data(ij)->data(iI);
          }
        }

        for (size_t ia = 0; ia != lena; ++ia) {
          for (auto& iter : cc_->det()->phib(k,l)) {
            size_t ibJ = iter.source;
            size_t ibI = iter.target;
            double sign = static_cast<double>(iter.sign);
            size_t iI = ibI + ia*lenb;
            size_t iJ = ibJ + ia*lenb;
            if ((iJ - offset) < dsize && iJ >= offset)
              e->element(iJ-offset, no) += sign * d->data(ij)->data(iI);
          }
        }

        if (i == l) {
          const int kj = k+j*norb_;
          for (size_t iJ = offset; iJ != offset+dsize; ++iJ) {
            e->element(iJ-offset, no) -= d->data(kj)->data(iJ);
          }
        }
        ++no;
      }
    }
  };

  // RDM3, RDM4 construction is multipassed and parallelized:
  //  (1) When ndet > 10000, (ndet < 10000 -> too small, almost no gain)
  //  and (2) When we have processes more than one
  //  OR  (3) When the number of words in <I|E_ij,kl|0> is larger than (10,10) case (635,040,000)
  const size_t ndet = cbra->det()->size();
  const size_t norb2 = norb_ * norb_;
  const size_t ijmax = 635040001 * 2;
  const size_t ijnum = ndet * norb2 * norb2;
  const size_t npass = ((mpi__->size() * 2 > ((ijnum-1)/ijmax + 1)) && (mpi__->size() != 1) && ndet > 10000) ? mpi__->size() * 2 : (ijnum-1) / ijmax + 1;
  const size_t nsize = (ndet-1) / npass + 1;
  Timer timer;
  if (npass > 1) {
    cout << "    * Third and Fock-weighted fourth order RDM (" << setw(2) << ist + 1 << "," << setw(2) << jst + 1 << ") evaluation" << endl;
    cout << "      will be done with " << npass << " passes" << endl;
  }

  rdm3->zero();
  rdm4f->zero();

  // multipassing through {I}
  for (size_t ipass = 0; ipass != npass; ++ipass) {
    if (ipass % mpi__->size() != mpi__->rank()) continue;

    const size_t ioffset = ipass * nsize;
    const size_t isize = (ipass != (npass - 1)) ? nsize : ndet - ioffset;
    const size_t halfsize = norb2 * (norb2 + 1) / 2;
    auto eket_half = make_shared<Matrix>(isize, halfsize, /*local=*/true);
    make_evec_half(dket, eket_half, isize, ioffset);
    shared_ptr<Matrix> ebra_half = eket_half;
    if (cbra != cket) {
      ebra_half = eket_half->clone();
      make_evec_half(dbra, ebra_half, isize, ioffset);
    }

    auto feket = make_shared<Matrix>(isize, norb2, /*local=*/true);
    feket->zero();
    {
      auto eket_full = make_shared<Matrix>(isize, norb2 * norb2, /*local=*/true);
      size_t no = 0;
      for (size_t kl = 0; kl != norb2; ++kl)
        for (size_t ij = kl; ij != norb2; ++ij) {
          size_t ijkl = ij + kl*norb2;
          size_t klij = kl + ij*norb2;
          for (size_t iI = 0; iI != isize; ++iI) {
            eket_full->element(iI, ijkl) = eket_half->element(iI, no);
            eket_full->element(iI, klij) = eket_half->element(iI, no);
          }
          ++no;
        }
      dgemv_("N", isize*norb_*norb_, norb_*norb_, 1.0, eket_full->data(), isize*norb_*norb_, fock->data(), 1, 0.0, feket->data(), 1);
    }

    auto dbram = make_shared<Matrix>(isize, norb2, /*local=*/true);
    for (size_t ij = 0; ij != norb2; ++ij)
      copy_n(&(dbra->data(ij)->data(ioffset)), isize, dbram->element_ptr(0, ij));

    // put in third-order RDM: <0|E_mn|I><I|E_ij < kl|0>
    auto tmp3 = make_shared<Matrix>(*dbram % *eket_half);
    auto tmp3_full = make_shared<Matrix>(norb2, norb2 * norb2, /*local=*/true);
    // <0|E_ij > kl|I>[I|mn|0]
    auto tmp4 = make_shared<Matrix>(*ebra_half % *feket);
    auto tmp4_full = make_shared<Matrix>(norb2 * norb2, norb2, /*local=*/true);

    for (size_t mn = 0; mn != norb2; ++mn) {
      size_t no = 0;
      for (size_t kl = 0; kl != norb2; ++kl) {
        for (size_t ij = kl; ij != norb2; ++ij) {
          size_t ijkl = ij + kl*norb2;
          size_t klij = kl + ij*norb2;
          tmp3_full->element(mn, ijkl) = tmp3->element(mn, no);
          tmp3_full->element(mn, klij) = tmp3->element(mn, no);
          tmp4_full->element(ijkl, mn) = tmp4->element(no, mn);
          tmp4_full->element(klij, mn) = tmp4->element(no, mn);
          ++no;
        }
      }
    }
    sort_indices<1,0,2,1,1,1,1>(tmp3_full->data(), rdm3->data(), norb_, norb_, norb2*norb2);
    sort_indices<1,0,3,2,4,1,1,1,1>(tmp4_full->data(), rdm4f->data(), norb_, norb_, norb_, norb_, norb2);
  }

  rdm3->allreduce();
  rdm4f->allreduce();

  if (npass > 1) {
    timer.tick_print("RDM evaluation (multipassing)");
  }

  {
    // then perform Eq. 49 of JCP 89 5803 (Werner's MRCI paper)
    // we assume that rdm2_[ist] is set
    for (int i0 = 0; i0 != norb_; ++i0)
      for (int i1 = 0; i1 != norb_; ++i1)
        for (int i2 = 0; i2 != norb_; ++i2)
          for (int i3 = 0; i3 != norb_; ++i3) {
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i2, i1, i0), norb_, rdm3->element_ptr(0, i3, i3, i2, i1, i0));
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i0, i3, i2), norb_, rdm3->element_ptr(0, i1, i3, i2, i1, i0));
          }
  }

  {
    // [0|E_ip,jl|0]
    auto frdm3 = make_shared<RDM<2>>(norb_);
    auto rdm3view = group(group(*rdm3, 4,6), 0,4);
    auto frdm3view = group(*frdm3, 0,4);
    contract(1.0, rdm3view, {0,1}, group(*fock, 0,2), {1}, 0.0, frdm3view, {0});

    // <0|E_ip,jn|0>f_nl (in this order)
    auto prdm2 = make_shared<RDM<2>>(norb_);
    auto prdm2v = group(*prdm2, 0,3);
    contract(1.0, group(*rdm2_->at(ist, jst), 0,3), {0,1}, *fock, {1,2}, 0.0, prdm2v, {0,2});

    // <0|E_ik,jl,on|0>f_np (in this order)
    auto prdm3 = make_shared<RDM<3>>(norb_);
    auto prdm3v = group(*prdm3, 0,5);
    contract(1.0, group(*rdm3, 0,5), {0,1}, *fock, {1,2}, 0.0, prdm3v, {0,2});

    sort_indices<2,0,1,1,1,-1,1>(prdm3->data(), rdm4f->data(), norb_*norb_, norb_*norb_, norb_*norb_);
    sort_indices<0,2,1,1,1,-1,1>(prdm3->data(), rdm4f->data(), norb_*norb_, norb_*norb_, norb_*norb_);

    auto prdm2t = prdm2->clone();
    sort_indices<1,0,0,1,1,1>(prdm2->data(), prdm2t->data(), norb_*norb_, norb_*norb_);
    for (int p = 0; p != norb_; ++p)
      for (int l = 0; l != norb_; ++l) {
        blas::ax_plus_y_n(-1.0, prdm2t->element_ptr(0,0,0,p), norb_*norb_*norb_, rdm4f->element_ptr(0,0,0,l,l,p));
        blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0,0,0,p),  norb_*norb_*norb_, rdm4f->element_ptr(0,0,0,l,l,p));
        for (int k = 0; k != norb_; ++k)
          for (int j = 0; j != norb_; ++j) {
            blas::ax_plus_y_n(-1.0, prdm2->element_ptr(0,p,j,l), norb_, rdm4f->element_ptr(0,k,j,l,k,p));
            blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0,p,j,l), norb_, rdm4f->element_ptr(0,k,j,l,k,p));
          }
      }
  }

  cc_->set_det(det_);

  return make_tuple(rdm3, rdm4f);
}
