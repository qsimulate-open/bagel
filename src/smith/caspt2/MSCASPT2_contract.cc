//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_contract.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_contract_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


shared_ptr<VectorB> MSCASPT2::MSCASPT2::contract_rdm_deriv_mat(shared_ptr<const CIWfn> ciwfn, int offset, int size, shared_ptr<const Matrix> fock,
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


shared_ptr<Queue> MSCASPT2::MSCASPT2::contract_rdm_deriv(shared_ptr<const CIWfn> ciwfn, shared_ptr<VectorB> bdata, int offset, int size, const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto contract = make_shared<Queue>();
  auto tensor900 = vector<shared_ptr<Tensor>>{deci};
  auto task900 = make_shared<Task900>(tensor900, reset);
  contract->add_task(task900);

#if 0
  auto tensor901 = vector<shared_ptr<Tensor>>{deci, rdm0deriv_, den0cit};
  auto task901 = make_shared<Task901>(tensor901, cindex);
  task901->add_dep(task900);
  contract->add_task(task901);

  auto tensor902 = vector<shared_ptr<Tensor>>{deci, rdm1deriv_, den1cit};
  auto task902 = make_shared<Task902>(tensor902, cindex);
  task902->add_dep(task900);
  contract->add_task(task902);

  auto tensor903 = vector<shared_ptr<Tensor>>{deci, rdm2deriv_, den2cit};
  auto task903 = make_shared<Task903>(tensor903, cindex);
  task903->add_dep(task900);
  contract->add_task(task903);
#endif

#if 0
  vector<IndexRange> I900_index = {ci_, active_, active_};
  auto I900 = make_shared<Tensor>(I900_index);
  auto tensor914 = vector<shared_ptr<Tensor>>{deci, I900};
  auto task914 = make_shared<Task914>(tensor914, cindex, ciwfn, bdata, offset, size);
  task914->add_dep(task900);
  contract->add_task(task914);

  vector<IndexRange> I901_index = {ci_, active_, active_};
  auto I901 = make_shared<Tensor>(I901_index);
#endif
  vector<IndexRange> I903_index = {active_, active_, active_, active_, active_, active_};
  auto I903 = make_shared<Tensor>(I903_index);
#if 0
  auto tensor915 = vector<shared_ptr<Tensor>>{I901, rdm2deriv_, den3cit, I903, deci};
  auto task915 = make_shared<Task915>(tensor915, cindex);
  task914->add_dep(task915);
  task915->add_dep(task900);
  contract->add_task(task915);
#endif

  auto tensor921 = vector<shared_ptr<Tensor>>{I903, den4cit, f1_};
  auto task921 = make_shared<Task921>(tensor921, cindex);
#if 0
  task915->add_dep(task921);
#endif
  task921->add_dep(task900);
  contract->add_task(task921);

  auto tensor916 = vector<shared_ptr<Tensor>>{deci, rdm2deriv_, den3cit, I903};
  auto task916 = make_shared<Task916>(tensor916, cindex);
  task916->add_dep(task921);
  task916->add_dep(task900);
  contract->add_task(task916);

#if 0
  vector<IndexRange> I902_index = {ci_, active_, active_};
  auto I902 = make_shared<Tensor>(I902_index);
  auto tensor917 = vector<shared_ptr<Tensor>>{I900, I901, I902, deci};
  auto task917 = make_shared<Task917>(tensor917, cindex);
  task917->add_dep(task900);
//  task917->add_dep(task916);
  task917->add_dep(task915);
  task914->add_dep(task917);
  contract->add_task(task917);

  auto tensor918 = vector<shared_ptr<Tensor>>{I902, rdm3fderiv_, den4cit, deci};
  auto task918 = make_shared<Task918>(tensor918, cindex);
  task917->add_dep(task918);
  task918->add_dep(task900);
  contract->add_task(task918);

#endif
#if 0
  auto tensor923 = vector<shared_ptr<Tensor>>{deci, rdm3fderiv_, den4cit};
  auto task923 = make_shared<Task923>(tensor923, cindex);
  task923->add_dep(task900);
  contract->add_task(task923);
#endif

  return contract;
}

#endif
