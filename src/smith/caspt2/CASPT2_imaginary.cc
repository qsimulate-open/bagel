//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_imaginary.cc
// Copyright (C) 2018 Toru Shiozaki
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
using namespace bagel::SMITH::CASPT2;


tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>,shared_ptr<RDM<3>>,shared_ptr<RDM<3>>> CASPT2::CASPT2::get_rdm() {
  const size_t nact  = info_->nact();

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  shared_ptr<RDM<3>> rdm3;
  shared_ptr<RDM<3>> rdm4f;

  // collect rdm1_
  {
    vector<IndexRange> o = rdm1_->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    auto d1 = make_shared<RDM<1>>(nact);
    for (auto& i1 : o[1].range())
      for (auto& i0 : o[0].range()) {
        auto input = rdm1_->get_block(i0, i1);
        for (size_t io1 = 0; io1 != i1.size(); ++io1)
          copy_n(&input[0 + i0.size() * io1], i0.size(), d1->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1));
      }
    rdm1 = d1->copy();
  }

  // collect rdm2_
  {
    vector<IndexRange>o = rdm2_->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    auto d2 = make_shared<RDM<2>>(nact);
    for (auto& i3 : o[3].range())
      for (auto& i2 : o[2].range())
        for (auto& i1 : o[1].range())
          for (auto& i0 : o[0].range()) {
            auto input = rdm2_->get_block(i0, i1, i2, i3);
            for (size_t io3 = 0; io3 != i3.size(); ++io3)
              for (size_t io2 = 0; io2 != i2.size(); ++io2)
                for (size_t io1 = 0; io1 != i1.size(); ++io1)
                  copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * io3))], i0.size(),
                         d2->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2, io3 + i3.offset() - off3));
          }
    rdm2 = d2->copy();
  }

  // collect rdm3_
  {
    vector<IndexRange>o = rdm3_->indexrange();
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
                auto input = rdm3_->get_block(i0, i1, i2, i3, i4, i5);
                for (size_t io5 = 0; io5 != i5.size(); ++io5)
                  for (size_t io4 = 0; io4 != i4.size(); ++io4)
                    for (size_t io3 = 0; io3 != i3.size(); ++io3)
                      for (size_t io2 = 0; io2 != i2.size(); ++io2)
                        for (size_t io1 = 0; io1 != i1.size(); ++io1)
                          copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * io5))))],
                                 i0.size(), d3->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                 io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5));
              }
    rdm3 = d3->copy();
  }

  // collect rdm4f_
  {
    vector<IndexRange>o = rdm4f_->indexrange();
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
                auto input = rdm4f_->get_block(i0, i1, i2, i3, i4, i5);
                for (size_t io5 = 0; io5 != i5.size(); ++io5)
                  for (size_t io4 = 0; io4 != i4.size(); ++io4)
                    for (size_t io3 = 0; io3 != i3.size(); ++io3)
                      for (size_t io2 = 0; io2 != i2.size(); ++io2)
                        for (size_t io1 = 0; io1 != i1.size(); ++io1)
                          copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * io5))))],
                                 i0.size(), d4->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                 io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5));
              }
    rdm4f = d4->copy();
  }

  return tie(rdm1, rdm2, rdm3, rdm4f);
}


void CASPT2::CASPT2::add_imaginary_shift(shared_ptr<Tensor> res, shared_ptr<Tensor> norm, const int istate) {
  const double shift2 = info_->shift() * info_->shift();
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  shared_ptr<RDM<3>> rdm3;
  shared_ptr<RDM<3>> rdm4f;

  tie(rdm1, rdm2, rdm3, rdm4f) = get_rdm();

  // a i b j case
  {
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            if (!norm->is_local(i0, i1, i2, i3) || !res->get_size(i0, i1, i2, i3)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i0, i1, i2, i3);
            unique_ptr<double[]> temp(new double[res->get_size(i0, i1, i2, i3)]);

            size_t iall = 0;
            for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
              for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1)
                for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
                  for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (eig_[j3] + eig_[j1] - eig_[j2] - eig_[j0]);
                  }

            res->add_block(temp, i0, i1, i2, i3);

          }
  }

  // a r b i case
  {
    for (auto& i0 : active_) {
      // trans is the transformation matrix
      assert(denom_->shalf_x());
      const size_t interm_size = denom_->shalf_x()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0) {
        unique_ptr<double[]> out(new double[I0.size()*interm_size]);
        for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
          copy_n(denom_->shalf_x()->element_ptr(0,j0-nclosed + i*nact), interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i0);

      for (auto& i3 : virt_) {
        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            const int blocksize = norm->get_size(i2, i3, i0, i1);
            unique_ptr<double[]> data0 = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> data2(new double[blocksize]);
            sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());

            // move to orthogonal basis
            unique_ptr<double[]> interm(new double[i1.size()*i2.size()*i3.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i2.size()*i3.size(), i0.size(),
                                        1.0, transp.get(), interm_size, data2.get(), i0.size(), 0.0, interm.get(), interm_size);

            size_t iall = 0;
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = 0; j0 != interm_size; ++j0, ++iall) {
                    interm[iall] *= shift2 / (denom_->denom_x(j0) + eig_[j3] + eig_[j1] - eig_[j2] - e0_);
                  }

            // move back to non-orthogonal basis
            unique_ptr<double[]> data3(new double[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size(), i1.size()*i2.size()*i3.size(), interm_size,
                                        1.0, transp.get(), interm_size, interm.get(), interm_size, 0.0, data3.get(), i0.size());
            unique_ptr<double[]> data4(new double[blocksize]);
            sort_indices<2,3,0,1,0,1,1,1>(data3, data4, i0.size(), i1.size(), i2.size(), i3.size());

            res->add_block(data4, i2, i3, i0, i1);
          }
        }
      }
    }
  }

  // a r b s case
  {
    for (auto& i2 : active_) {
    for (auto& i0 : active_) {
      // trans is the transformation matrix
      assert(denom_->shalf_xx());
      const size_t interm_size = denom_->shalf_xx()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I2) {
        unique_ptr<double[]> out(new double[I0.size()*I2.size()*interm_size]);
        for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
          for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
            copy_n(denom_->shalf_xx()->element_ptr(0,(j0-nclosed)+(j2-nclosed)*nact + i*nact*nact),
                   interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i0, i2);

      for (auto& i3 : virt_) {
        for (auto& i1 : virt_) {
          if (!norm->is_local(i0, i1, i2, i3) || !res->get_size(i0, i1, i2, i3)) continue;
          const size_t blocksize = norm->get_size(i0, i1, i2, i3);
          unique_ptr<double[]> data0 = norm->get_block(i0, i1, i2, i3);
          unique_ptr<double[]> data1(new double[blocksize]);
          // sort. Active indices run faster
          sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i1.size()*i3.size()*interm_size]);

          // move to orthogonal basis
          btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i3.size(), i0.size()*i2.size(),
                                      1.0, transp.get(), interm_size, data1.get(), i0.size()*i2.size(), 0.0, interm.get(), interm_size);

          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j02 = 0; j02 != interm_size; ++j02, ++iall)
                interm[iall] *= shift2 / (denom_->denom_xx(j02) + eig_[j3] + eig_[j1] - e0_);

          btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), interm_size,
                                      0.5, transp.get(), interm_size, interm.get(), interm_size, 0.0, data1.get(), i0.size()*i2.size());

          // sort back to the original order
          unique_ptr<double[]> data2(new double[blocksize]);
          sort_indices<0,2,1,3,0,1,1,1>(data1, data2, i0.size(), i2.size(), i1.size(), i3.size());
          res->add_block(data2, i0, i1, i2, i3);
        }
      }
    }
    }
  }
#if 0
#if 0
  // a r b i case
  {
    auto dx = make_shared<RDM<1>>(nact);
    auto workv = btas::group(*dx, 0,2);
    auto rdm2v = btas::group(btas::group(*rdm2,2,4),0,2);
    btas::contract(1.0, rdm2v, {0,1}, btas::group(*fockact_,0,2), {1}, 0.0, workv, {0});
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> temp(new double[res->get_size(i2, i3, i0, i1)]);

            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
              for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3)
                for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
                  for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (dx->element(j0-nclosed,j0-nclosed) - eig_[j2] + eig_[j3] + eig_[j1] - e0_);
                  }
            res->add_block(temp, i2, i3, i0, i1);
          }
  }
#endif
  // a r b s case
  {
    // this term is either wrong or I am misunderstanding something
    auto dxx = make_shared<RDM<2>>(nact);
    auto workv = btas::group(*dxx, 0,4);
    auto rdm3v = btas::group(btas::group(*rdm3,4,6),0,4);
    btas::contract(1.0, rdm3v, {0,1}, btas::group(*fockact_,0,2), {1}, 0.0, workv, {0});
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i3 : virt_)
          for (auto& i1 : virt_) {
            if (!norm->is_local(i0, i1, i2, i3) || !res->get_size(i0, i1, i2, i3)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i0, i1, i2, i3);
            unique_ptr<double[]> temp(new double[res->get_size(i0, i1, i2, i3)]);

            size_t iall = 0;
            for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
              for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1)
                for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
                  for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (dxx->element(j0-nclosed,j2-nclosed,j0-nclosed,j2-nclosed) + eig_[j3] + eig_[j1] - e0_ );
                  }
            res->add_block(temp, i0, i1, i2, i3);
          }
  }

  // a i r j case
  // responsible for different sth
  {
    const size_t dim = nact;
    auto dh = make_shared<Matrix>(dim, dim);
    {
      auto D2h = make_shared<RDM<2>>(*rdm2);
      D2h->scale(-1.0);
      for (int i0 = 0; i0 != nact; ++i0)
        for (int i1 = 0; i1 != nact; ++i1)
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i3 = 0; i3 != nact; ++i3) {
              if (i1 == i2 && i0 == i3) D2h->element(i0, i1, i2, i3) += 2.0;
              if (i1 == i2) D2h->element(i0, i1, i2, i3) -= rdm1->element(i0, i3);
              if (i0 == i3) D2h->element(i0, i1, i2, i3) -= rdm1->element(i1, i2);
              if (i0 == i1) D2h->element(i0, i1, i2, i3) += 2.0 * rdm1->element(i2, i3);
            }
      auto workv = btas::group(*dh, 0,2);
      auto rdm2v = btas::group(btas::group(*D2h,2,4),0,2);
      btas::contract(1.0, rdm2v, {0,1}, btas::group(*fockact_,0,2), {1}, 0.0, workv, {0});
    }

    for (auto& i3 : active_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> temp(new double[res->get_size(i2, i3, i0, i1)]);

            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
              for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3)
                for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
                  for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (dh->element(j3-nclosed,j3-nclosed) - eig_[j2] + eig_[j1] - eig_[j0] - e0_);
                  }

            res->add_block(temp, i2, i3, i0, i1);
          }
  }

  // r i s j case
  // seems right
  {
    const size_t dim = nact * nact;
    auto dhh = make_shared<Matrix>(dim, dim);
    {
      auto D3hh = make_shared<RDM<3>>(*rdm3);
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          for (int k = 0; k != nact; ++k)
            for (int l = 0; l != nact; ++l)
              for (int m = 0; m != nact; ++m)
                for (int n = 0; n != nact; ++n) {
                  if (l == n) D3hh->element(i,j,k,l,m,n) += rdm2->element(j,m,i,k);
                  if (l == m) D3hh->element(i,j,k,l,m,n) += rdm2->element(i,n,j,k);
                  if (i == n) D3hh->element(i,j,k,l,m,n) += rdm2->element(i,m,l,k);
                  if (j == n && l == m) D3hh->element(i,j,k,l,m,n) += rdm1->element(i,k);
                  if (j == m) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm2->element(i,n,l,k);
                  if (j == m && l == n) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm1->element(i,k);
                  if (j == k) D3hh->element(i,j,k,l,m,n) += rdm2->element(i,n,l,m);
                  if (j == k && l == n) D3hh->element(i,j,k,l,m,n) += rdm1->element(i,m);
                  if (j == k && l == m) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm1->element(i,n);
                  if (i == n) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm2->element(j,m,l,k);
                  if (i == n && l == m) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm1->element(j,k);
                  if (j == k && i == n) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm1->element(l,m);
                  if (j == k && i == n && l == m) D3hh->element(i,j,k,l,m,n) += 4.0;
                  if (i == m) D3hh->element(i,j,k,l,m,n) += rdm2->element(i,n,l,k);
                  if (l == n && i == m) D3hh->element(i,j,k,l,m,n) += rdm1->element(j,k);
                  if (j == k && i == m) D3hh->element(i,j,k,l,m,n) += rdm1->element(l,n);
                  if (j == k && i == m && l == n) D3hh->element(i,j,k,l,m,n) -= 2.0;
                  if (i == k) D3hh->element(i,j,k,l,m,n) += rdm2->element(l,n,j,m);
                  if (i == k && l == n) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm1->element(j,m);
                  if (i == k && l == m) D3hh->element(i,j,k,l,m,n) += rdm1->element(j,n);
                  if (j == n && i == k) D3hh->element(i,j,k,l,m,n) += rdm1->element(i,m);
                  if (i == k && j == n && l == m) D3hh->element(i,j,k,l,m,n) -= 2.0;
                  if (i == k && j == m) D3hh->element(i,j,k,l,m,n) += rdm1->element(l,n);
                  if (j == m && l == n && i == k) D3hh->element(i,j,k,l,m,n) += 4.0;
                  if (i == n && j == m) D3hh->element(i,j,k,l,m,n) += 4.0 * rdm1->element(l,k);
                  if (i == m && j == n) D3hh->element(i,j,k,l,m,n) -= 2.0 * rdm1->element(l,k);
                }
      auto workv = btas::group(*dhh, 0,2);
      auto rdm3v = btas::group(btas::group(*D3hh,4,6),0,4);
      btas::contract(1.0, rdm3v, {0,1}, btas::group(*fockact_,0,2), {1}, 0.0, workv, {0});
    }
    for (auto& i3 : active_)
      for (auto& i1 : active_)
        for (auto& i2 : closed_)
          for (auto& i0 : closed_) {
            if (!norm->is_local(i0, i1, i2, i3) || !res->get_size(i0, i1, i2, i3)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i0, i1, i2, i3);
            unique_ptr<double[]> temp(new double[res->get_size(i0, i1, i2, i3)]);

            size_t iall = 0;
            for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
              for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1)
                for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
                  for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (dhh->element(j1-nclosed+(j3-nclosed)*nact, j1-nclosed+(j3-nclosed)*nact) - eig_[j2] - eig_[j0] - e0_);
                  }

            res->add_block(temp, i0, i1, i2, i3);
          }
  }

  // a i r s case
  // wrong from here
  {
    const size_t dim = nact * nact;
    auto dxh2 = make_shared<Matrix>(dim, dim);
    {
      auto D3xh2 = make_shared<RDM<3>>(*rdm3);
      D3xh2->scale(-1.0);
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          for (int k = 0; k != nact; ++k)
            for (int l = 0; l != nact; ++l)
              for (int m = 0; m != nact; ++m)
                for (int n = 0; n != nact; ++n) {
                  if (n == i) D3xh2->element(i,j,k,l,m,n) -= rdm2->element(m,j,k,l);
                  if (j == m) D3xh2->element(i,j,k,l,m,n) -= rdm2->element(i,n,k,l);
                  if (n == i && j == m) D3xh2->element(i,j,k,l,m,n) += 2.0 * rdm1->element(k,l);
                  if (i == j) D3xh2->element(i,j,k,l,m,n) += 2.0 * rdm2->element(m,n,k,l);
                }
      auto workv = btas::group(*dxh2, 0,2);
      auto rdm3v = btas::group(btas::group(*D3xh2,4,6),0,4);
      btas::contract(1.0, rdm3v, {0,1}, btas::group(*fockact_,0,2), {1}, 0.0, workv, {0});
    }
    for (auto& i0 : closed_)
      for (auto& i1 : virt_)
        for (auto& i2 : active_)
          for (auto& i3 : active_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> temp(new double[res->get_size(i2, i3, i0, i1)]);

            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
              for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3)
                for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
                  for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (dxh2->element(j2-nclosed+(j2-nclosed)*nact, j3-nclosed+(j3-nclosed)*nact) + eig_[j1] - eig_[j0] - e0_);
//                    cout << " " << j0 << " " << j1 << " " << j2 << " " << j3 << "  = " << dxh1->element(j2+j2*nact, j3+j3*nact) + eig_[j1] - eig_[j0] - e0_ << endl;
                  }

            res->add_block(temp, i2, i3, i0, i1);
          }
  }

//  cout << " a r s i case " << endl;
  // a r s i case
  {
    const size_t dim = nact * nact;
    auto dxh1 = make_shared<Matrix>(dim, dim);
    {
      auto D3xh1 = make_shared<RDM<3>>(*rdm3);
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          for (int k = 0; k != nact; ++k)
            for (int l = 0; l != nact; ++l)
              for (int m = 0; m != nact; ++m)
                for (int n = 0; n != nact; ++n) {
                  if (n == k) D3xh1->element(i,j,k,l,m,n) += rdm2->element(i,j,m,l);
                  if (j == m) D3xh1->element(i,j,k,l,m,n) += rdm2->element(i,n,k,l);
                  if (j == m && n == k) D3xh1->element(i,j,k,l,m,n) += rdm1->element(i,l);
                  if (j == k) D3xh1->element(i,j,k,l,m,n) += rdm2->element(m,n,i,l);
                }
      auto workv = btas::group(*dxh1, 0,2);
      auto rdm3v = btas::group(btas::group(*D3xh1,4,6),0,4);
      btas::contract(1.0, rdm3v, {0,1}, btas::group(*fockact_,0,2), {1}, 0.0, workv, {0});
    }
    for (auto& i0 : closed_)
      for (auto& i1 : virt_)
        for (auto& i2 : active_)
          for (auto& i3 : active_) {
            if (!norm->is_local(i0, i3, i2, i1) || !res->get_size(i0, i3, i2, i1)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i0, i3, i2, i1);
            unique_ptr<double[]> temp(new double[res->get_size(i0, i3, i2, i1)]);

            size_t iall = 0;
            for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
              for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1, ++iall) {
                    temp[iall] = ncurrent[iall] * shift2 / (dxh1->element(j2-nclosed+(j2-nclosed)*nact, j3-nclosed+(j3-nclosed)*nact) + eig_[j1] - eig_[j0] - e0_);
//                    cout << " " << j0 << " " << j1 << " " << j2 << " " << j3 << "  = " << dxh2->element(j2+j2*nact, j3+j3*nact) + eig_[j1] - eig_[j0] - e0_ << endl;
                  }

            res->add_block(temp, i0, i3, i2, i1);
          }
  }

  // a r s t case
  {
    auto dxxh = make_shared<RDM<3>>(*rdm4f);
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        for (int k = 0; k != nact; ++k)
          for (int l = 0; l != nact; ++l)
            for (int m = 0; m != nact; ++m)
              for (int n = 0; n != nact; ++n) {
                dxxh->element(i,j,k,l,m,n) += rdm2->element(i,l,m,n) * fockact_->element(j,k);
                for (int o = 0; o != nact; ++o) {
                  dxxh->element(i,j,k,l,m,n) += rdm3->element(i,j,o,l,m,n) * fockact_->element(o,k);
                  dxxh->element(i,j,k,l,m,n) += rdm3->element(i,o,k,l,m,n) * fockact_->element(j,o);
                  for (int p = 0; p != nact; ++p)
                    if (j == k) dxxh->element(i,j,k,l,m,n) += rdm3->element(o,p,i,l,m,n) * fockact_->element(o,p);
                }
              }
    for (auto& i0 : active_)
      for (auto& i1 : virt_)
        for (auto& i2 : active_)
          for (auto& i3 : active_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> temp(new double[res->get_size(i2, i3, i0, i1)]);

            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
              for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3)
                for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
                  for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1, ++iall) {
//                    temp[iall] = ncurrent[iall] * shift2 / (dxxh->element(j0-nclosed,j0-nclosed,j2-nclosed,j2-nclosed,j3-nclosed,j3-nclosed) + eig_[j1] - e0_);
                    temp[iall] = ncurrent[iall] * shift2 / (dxxh->element(j0-nclosed,j2-nclosed,j3-nclosed,j0-nclosed,j2-nclosed,j3-nclosed) + eig_[j1] - e0_);
                  }

            res->add_block(temp, i2, i3, i0, i1);
          }
  }
 
  // r i s t case
  {
    auto dxhh = make_shared<RDM<3>>(*rdm4f);
    dxhh->scale(-1.0);
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        for (int k = 0; k != nact; ++k)
          for (int l = 0; l != nact; ++l)
            for (int m = 0; m != nact; ++m)
              for (int n = 0; n != nact; ++n) {
                dxhh->element(i,j,k,l,m,n) -= rdm2->element(i,j,k,m) * fockact_->element(l,m);
                dxhh->element(i,j,k,l,m,n) += 2.0 * rdm2->element(i,j,m,n) * fockact_->element(l,k);
                if (j == m) dxhh->element(i,j,k,l,m,n) += 2.0 * rdm1->element(i,n) * fockact_->element(l,k);
                if (j == k) dxhh->element(i,j,k,l,m,n) -= rdm1->element(i,n) * fockact_->element(l,m);
                dxhh->element(i,j,k,l,m,n) -= rdm2->element(k,l,i,n) * fockact_->element(j,m);
                dxhh->element(i,j,k,l,m,n) -= rdm2->element(i,l,m,n) * fockact_->element(j,k);
                if (l == m) dxhh->element(i,j,k,l,m,n) -= rdm1->element(i,n) * fockact_->element(j,k);
                if (l == k && j == m) dxhh->element(i,j,k,l,m,n) += 2.0 * rdm1->element(i,n);
                for (int o = 0; o != nact; ++o) {
                  dxhh->element(i,j,k,l,m,n) -= rdm3->element(o,n,i,j,k,l) * fockact_->element(o,m);
                  dxhh->element(i,j,k,l,m,n) -= rdm3->element(o,l,m,n,i,j) * fockact_->element(o,k);
                  if (l == m) dxhh->element(i,j,k,l,m,n) -= rdm2->element(o,n,i,j) * fockact_->element(o,k);
                  if (l == k) dxhh->element(i,j,k,l,m,n) += 2.0 * rdm2->element(o,n,i,j) * fockact_->element(o,m);
                  dxhh->element(i,j,k,l,m,n) -= rdm3->element(m,n,i,j,k,o) * fockact_->element(l,o);
                  if (j == m) dxhh->element(i,j,k,l,m,n) -= rdm2->element(o,l,i,n) * fockact_->element(o,k);
                  if (j == m) dxhh->element(i,j,k,l,m,n) -= rdm2->element(i,n,k,o) * fockact_->element(l,o);
                  if (j == k) dxhh->element(i,j,k,l,m,n) -= rdm2->element(o,n,i,l) * fockact_->element(o,m);
                  if (j == k) dxhh->element(i,j,k,l,m,n) -= rdm2->element(m,n,i,o) * fockact_->element(l,o);
                  dxhh->element(i,j,k,l,m,n) -= rdm3->element(k,l,m,n,i,o) * fockact_->element(j,o);
                  if (l == m) dxhh->element(i,j,k,l,m,n) -= rdm2->element(k,n,i,o) * fockact_->element(j,o);
                  if (l == k) dxhh->element(i,j,k,l,m,n) += 2.0 * rdm2->element(m,n,i,o) * fockact_->element(j,o);
                  for (int p = 0; p != nact; ++p) {
                    if (l == m) dxhh->element(i,j,k,l,m,n) -= rdm3->element(k,n,i,j,o,p) * fockact_->element(o,p);
                    if (k == l) dxhh->element(i,j,k,l,m,n) += rdm3->element(m,n,i,j,o,p) * fockact_->element(o,p);
                    if (j == m) dxhh->element(i,j,k,l,m,n) -= rdm3->element(k,l,i,n,o,p) * fockact_->element(o,p);
                    if (j == k) dxhh->element(i,j,k,l,m,n) -= rdm3->element(m,n,i,l,o,p) * fockact_->element(o,p);
                    if (j == m && l == k) dxhh->element(i,j,k,l,m,n) += 2.0 * rdm2->element(i,n,o,p) * fockact_->element(o,p);
                    if (j == k && l == m) dxhh->element(i,j,k,l,m,n) -= rdm2->element(i,n,o,p) * fockact_->element(o,p);
                  }
                }
              }

    for (auto& i3 : active_)
      for (auto& i1 : active_)
        for (auto& i2 : closed_)
          for (auto& i0 : active_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            unique_ptr<double[]> ncurrent = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> temp(new double[res->get_size(i2, i3, i0, i1)]);

            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset() + i2.size(); ++j2)
              for (int j3 = i3.offset(); j3 != i3.offset() + i3.size(); ++j3)
                for (int j0 = i0.offset(); j0 != i0.offset() + i0.size(); ++j0)
                  for (int j1 = i1.offset(); j1 != i1.offset() + i1.size(); ++j1, ++iall) {
//                    temp[iall] = ncurrent[iall] * shift2 / (dxhh->element(j0-nclosed,j0-nclosed,j1-nclosed,j1-nclosed,j3-nclosed,j3-nclosed) - eig_[j2] - e0_);
                    temp[iall] = ncurrent[iall] * shift2 / (dxhh->element(j0-nclosed,j1-nclosed,j3-nclosed,j0-nclosed,j1-nclosed,j3-nclosed) - eig_[j2] - e0_);
                  }

            res->add_block(temp, i2, i3, i0, i1);
          }
  }
#endif

}

#endif
