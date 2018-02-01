//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_contract.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/caspt2/CASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


tuple<shared_ptr<double>,shared_ptr<RDM<1>>,shared_ptr<RDM<2>>,shared_ptr<RDM<3>>,shared_ptr<RDM<3>>> CASPT2::CASPT2::feed_denci() {
  const size_t nact  = info_->nact();

  shared_ptr<double> den0cirdm;
  shared_ptr<RDM<1>> den1cirdm;
  shared_ptr<RDM<2>> den2cirdm;
  shared_ptr<RDM<3>> den3cirdm;
  shared_ptr<RDM<3>> den4cirdm;

  den0cit = den0ci;
  den1cit = den1ci;
  den2cit = den2ci;
  den3cit = den3ci;
  den4cit = den4ci;
  mpi__->barrier();

  // collect den0ci
  {
    unique_ptr<double[]> d0data = den0cit->get_block();
    auto d0 = make_shared<double>(d0data[0]);
    den0cirdm = d0;
  }

  // collect den1ci
  {
    vector<IndexRange> o = den1cit->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    auto d1 = make_shared<RDM<1>>(nact);
    for (auto& i1 : o[1].range())
      for (auto& i0 : o[0].range()) {
        auto input = den1cit->get_block(i0, i1);
        for (size_t io1 = 0; io1 != i1.size(); ++io1)
          copy_n(&input[0+i0.size()*io1], i0.size(), d1->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1));
      }
    den1cirdm = d1->copy();
  }

  // collect den2ci
  {
    vector<IndexRange>o = den2cit->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    auto d2 = make_shared<RDM<2>>(nact);
    for (auto& i3 : o[3].range())
      for (auto& i2 : o[2].range())
        for (auto& i1 : o[1].range())
          for (auto& i0 : o[0].range()) {
            auto input = den2cit->get_block(i0, i1, i2, i3);
            for (size_t io3 = 0; io3 != i3.size(); ++io3)
              for (size_t io2 = 0; io2 != i2.size(); ++io2)
                for (size_t io1 = 0; io1 != i1.size(); ++io1)
                  copy_n(&input[0+i0.size()*(io1+i1.size()*(io2+i2.size()*io3))], i0.size(), d2->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2, io3 + i3.offset() - off3));
          }
    den2cirdm = d2->copy();
  }

  // collect den3ci
  {
    vector<IndexRange>o = den3cit->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    const int off4 = o[4].front().offset();
    const int off5 = o[5].front().offset();
    auto d3 = make_shared<RDM<3>>(nact);
    for (auto& i5 : o[5].range())
      for (auto& i4 : o[4].range())
        for (auto& i3 : o[3].range())
          for (auto& i2 : o[2].range())
            for (auto& i1 : o[1].range())
              for (auto& i0 : o[0].range()) {
                auto input = den3cit->get_block(i0, i1, i2, i3, i4, i5);
                for (size_t io5 = 0; io5 != i5.size(); ++io5)
                  for (size_t io4 = 0; io4 != i4.size(); ++io4)
                    for (size_t io3 = 0; io3 != i3.size(); ++io3)
                      for (size_t io2 = 0; io2 != i2.size(); ++io2)
                        for (size_t io1 = 0; io1 != i1.size(); ++io1)
                          copy_n(&input[0+i0.size()*(io1+i1.size()*(io2+i2.size()*(io3+i3.size()*(io4+i4.size()*io5))))],
                              i0.size(), d3->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2, io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5));
              }
    den3cirdm = d3->copy();
  }

  // collect den4ci
  {
    vector<IndexRange>o = den4cit->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    const int off4 = o[4].front().offset();
    const int off5 = o[5].front().offset();
    auto d4 = make_shared<RDM<3>>(nact);
    for (auto& i5 : o[5].range())
      for (auto& i4 : o[4].range())
        for (auto& i3 : o[3].range())
          for (auto& i2 : o[2].range())
            for (auto& i1 : o[1].range())
              for (auto& i0 : o[0].range()) {
                auto input = den4cit->get_block(i0, i1, i2, i3, i4, i5);
                for (size_t io5 = 0; io5 != i5.size(); ++io5)
                  for (size_t io4 = 0; io4 != i4.size(); ++io4)
                    for (size_t io3 = 0; io3 != i3.size(); ++io3)
                      for (size_t io2 = 0; io2 != i2.size(); ++io2)
                        for (size_t io1 = 0; io1 != i1.size(); ++io1)
                          copy_n(&input[0+i0.size()*(io1+i1.size()*(io2+i2.size()*(io3+i3.size()*(io4+i4.size()*io5))))],
                              i0.size(), d4->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2, io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5));
              }
    den4cirdm = d4->copy();
  }

  return tie(den0cirdm, den1cirdm, den2cirdm, den3cirdm, den4cirdm);
}


shared_ptr<VectorB> CASPT2::CASPT2::contract_rdm_deriv(shared_ptr<const CIWfn> ciwfn, int offset, int size, shared_ptr<const Matrix> fock,
    shared_ptr<const VectorB> rdm0deriv, shared_ptr<const Matrix> rdm1deriv, shared_ptr<const Matrix> rdm2deriv, shared_ptr<const Matrix> rdm3deriv,
    shared_ptr<const double> den0cirdmt, shared_ptr<const RDM<1>> den1cirdmt, shared_ptr<const RDM<2>> den2cirdmt, shared_ptr<const RDM<3>> den3cirdmt, shared_ptr<const RDM<3>> den4cirdmt) {
  const size_t ndet = ci_deriv_->data(0)->size();
  const size_t nact  = info_->nact();
  auto out = make_shared<VectorB>(ndet);

  // rdm0deriv contraction
  {
    blas::ax_plus_y_n(den0cirdmt.get()[0], rdm0deriv->data(), size, out->data()+offset);
  }

  // rdm1deriv contraction
  {
    const size_t nact2 = nact * nact;
    dgemv_("N", size, nact2, 1.0, rdm1deriv->data(), size, den1cirdmt->data(), 1, 1.0, out->data()+offset, 1);
  }

  // rdm2deriv contraction
  {
    const size_t nact2 = nact * nact;
    const size_t nact4 = nact2 * nact2;
    dgemv_("N", size, nact4, 1.0, rdm2deriv->data(), size, den2cirdmt->data(), 1, 1.0, out->data()+offset, 1);
  }

  // rdm3, 4 deriv contraction
  {
    const size_t nact2 = nact * nact;
    const size_t nact4 = nact2 * nact2;
    const size_t nact5 = nact4 * nact;
    shared_ptr<RDM<3>> den3cif = den4cirdmt->clone();
    dgemm_("N", "T", nact5, nact, nact, -1.0, den4cirdmt->data(), nact5, fock->data(), nact, 1.0, den3cif->data(), nact5);

    const size_t nact6 = nact4 * nact2;
    blas::ax_plus_y_n(1.0, den3cirdmt->data(), nact6, den3cif->data());

    {
      shared_ptr<RDM<2>> den3cis = den2cirdmt->clone();
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          for (int l = 0; l != nact; ++l)
            for (int m = 0; m != nact; ++m)
              for (int n = 0; n != nact; ++n) {
                den3cis->element(i,l,m,n) += den3cif->element(i,j,j,l,m,n);
                den3cis->element(i,l,m,n) += den3cif->element(m,j,i,l,j,n);
              }
      dgemv_("N", size, nact4, -1.0, rdm2deriv->data(), size, den3cis->data(), 1, 1.0, out->data()+offset, 1);
    }

    {
      shared_ptr<RDM<2>> den3cis = den2cirdmt->clone();
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          for (int l = 0; l != nact; ++l)
            for (int m = 0; m != nact; ++m)
              for (int n = 0; n != nact; ++n) {
                den3cis->element(i,l,m,n) += den4cirdmt->element(j,n,i,l,m,j);
                den3cis->element(i,l,m,n) += den4cirdmt->element(i,l,j,n,m,j);
              }
      dgemv_("N", size, nact4, -1.0, rdm3deriv->data(), size, den3cis->data(), 1, 1.0, out->data()+offset, 1);
    }

    {
      shared_ptr<Matrix> dtensor = rdm1deriv->clone();
      // TODO these operations should be made half (nact4 -> nunique)
      dgemm_("N", "T", size, nact2, nact4, 1.0, rdm2deriv->data(), size, den3cif->data(), nact2, 1.0, dtensor->data(), size);
      dgemm_("N", "N", size, nact2, nact4, 1.0, rdm3deriv->data(), size, den4cirdmt->data(), nact4, 1.0, dtensor->data(), size);

      {
        const size_t lena = ciwfn->det()->lena();
        const size_t lenb = ciwfn->det()->lenb();

        for (size_t ij = 0; ij != nact2; ++ij) {
          const size_t j = ij/nact;
          const size_t i = ij-j*nact;

          for (auto& iter : ciwfn->det()->phia(i,j)) {
            size_t iaJ = iter.source;
            size_t iaK = iter.target;
            double sign = static_cast<double>(iter.sign);
            for (size_t ib = 0; ib != lenb; ++ib) {
              size_t iK = ib+iaK*lenb;
              size_t iJ = ib+iaJ*lenb;
              if ((iK - offset) < size && iK >= offset)
                (*out)[iJ] += sign * dtensor->element(iK-offset,ij);
            }
          }

          for (size_t ia = 0; ia != lena; ++ia) {
            for (auto& iter : ciwfn->det()->phib(i,j)) {
              size_t ibJ = iter.source;
              size_t ibK = iter.target;
              double sign = static_cast<double>(iter.sign);
              size_t iK = ibK+ia*lenb;
              size_t iJ = ibJ+ia*lenb;
              if ((iK - offset) < size && iK >= offset)
                (*out)[iJ] += sign * dtensor->element(iK-offset,ij);
            }
          }
        }
      }
    }
  }

  return out;
}


#endif
