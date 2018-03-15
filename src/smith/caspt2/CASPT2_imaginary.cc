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


void CASPT2::CASPT2::add_imaginary_shift(shared_ptr<Tensor> res, shared_ptr<const Tensor> norm, const int istate) {
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

  // a r b i case
  {
    for (auto& i0 : active_) {
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
                    interm[iall] *= (shift2 / (denom_->denom_x(j0) + eig_[j3] + eig_[j1] - eig_[j2] - e0_));
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
#if 1

  // a i r j case
  {
    for (auto& i3 : active_) {
      // trans is the transformation matrix
      assert(denom_->shalf_h());
      const size_t interm_size = denom_->shalf_x()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I3) {
        unique_ptr<double[]> out(new double[I3.size()*interm_size]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
          copy_n(denom_->shalf_h()->element_ptr(0,j3-nclosed + i*nact), interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i3);

      for (auto& i2 : closed_) {
        for (auto& i1 : virt_) {
          for (auto& i0 : closed_) {
            if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
            const size_t blocksize = norm->get_size(i2, i3, i0, i1);

            unique_ptr<double[]> data0 = norm->get_block(i2, i3, i0, i1);
            unique_ptr<double[]> data2(new double[blocksize]);
            sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());

            // move to orthogonal basis
            unique_ptr<double[]> interm(new double[i0.size()*i1.size()*i2.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size()*i2.size(), interm_size, i3.size(),
                                        1.0, data2.get(), i0.size()*i1.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size()*i2.size());

            size_t iall = 0;
            for (int j3 = 0; j3 != interm_size; ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] *= shift2 / (denom_->denom_h(j3) - eig_[j2] + eig_[j1] - eig_[j0] - e0_);

            // move back to non-orthogonal basis
            unique_ptr<double[]> data3(new double[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size()*i2.size(), i3.size(), interm_size,
                                        1.0, interm.get(), i0.size()*i1.size()*i2.size(), transp.get(), interm_size, 0.0, data3.get(), i0.size()*i1.size()*i2.size());

            unique_ptr<double[]> data4(new double[blocksize]);
            sort_indices<2,3,0,1,0,1,1,1>(data3, data4, i0.size(), i1.size(), i2.size(), i3.size());

            res->add_block(data4, i2, i3, i0, i1);
          }
        }
      }
    }
  }

  // r i s j case
  {
    for (auto& i3 : active_) {
    for (auto& i1 : active_) {
      assert(denom_->shalf_hh());
      const size_t interm_size = denom_->shalf_hh()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I1, const Index& I3) {
        unique_ptr<double[]> out(new double[I1.size()*I3.size()*interm_size]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
          for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
            copy_n(denom_->shalf_hh()->element_ptr(0,(j1-nclosed)+(j3-nclosed)*nact + i*nact*nact),
                   interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i1, i3);

      for (auto& i2 : closed_) {
        for (auto& i0 : closed_) {
          if (!norm->is_local(i0, i1, i2, i3) || !res->get_size(i0, i1, i2, i3)) continue;
          const size_t blocksize = norm->get_size(i0, i1, i2, i3);
          // data0 is the source area
          unique_ptr<double[]> data0 = norm->get_block(i0, i1, i2, i3);
          unique_ptr<double[]> data1(new double[blocksize]);
          // sort. Active indices run slower
          sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i0.size()*i2.size()*interm_size]);

          // move to orthogonal basis
          btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i2.size(), interm_size, i1.size()*i3.size(),
                                      1.0, data1.get(), i0.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i2.size());

          size_t iall = 0;
          for (int j13 = 0; j13 != interm_size; ++j13)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                interm[iall] *= shift2 / (denom_->denom_hh(j13) - eig_[j2] - eig_[j0] - e0_);

          // move back to non-orthogonal basis
          // factor of 0.5 due to the factor in the overlap
          btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), interm_size,
                                      0.5, interm.get(), i0.size()*i2.size(), transp.get(), interm_size, 0.0, data1.get(), i0.size()*i2.size());

          // sort back to the original order
          unique_ptr<double[]> data2(new double[blocksize]);
          sort_indices<0,2,1,3,0,1,1,1>(data1, data2, i0.size(), i2.size(), i1.size(), i3.size());
          res->add_block(data2, i0, i1, i2, i3);
        }
      }
    }
    }
  }

  // a i r s & a r s i case
  {
    for (auto& i3 : active_) {
    for (auto& i2 : active_) {
      assert(denom_->shalf_xh());
      const size_t interm_size = denom_->shalf_xh()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I2, const Index& I3) {
        unique_ptr<double[]> out(new double[I2.size()*I3.size()*interm_size*2]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
          for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
            copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclosed)+(j3-nclosed)*nact + 2*i*nact*nact),
                   interm_size, out.get()+interm_size*k);
            copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclosed)+(j3-nclosed)*nact + (2*i+1)*nact*nact),
                   interm_size, out.get()+interm_size*(k+I2.size()*I3.size()));
          }
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i2, i3);

      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
          const size_t blocksize = norm->get_size(i2, i3, i0, i1);
          unique_ptr<double[]> data0 = norm->get_block(i2, i3, i0, i1);
          unique_ptr<double[]> data1 = norm->get_block(i0, i3, i2, i1);

          unique_ptr<double[]> data2(new double[blocksize*2]);
          // sort. Active indices run slower
          sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
          sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i0.size()*i1.size()*interm_size]);

          // move to orthogonal basis
          btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size(), interm_size, i2.size()*i3.size()*2,
                                      1.0, data2.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size());

          size_t iall = 0;
          for (int j23 = 0; j23 != interm_size; ++j23)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                interm[iall] *= shift2 / (denom_->denom_xh(j23) + eig_[j1] - eig_[j0] - e0_);

          // move back to non-orthogonal basis
          btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2.size()*i3.size()*2, interm_size,
                                      1.0, interm.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size());

          // sort back to the original order
          unique_ptr<double[]> data3(new double[blocksize]);
          unique_ptr<double[]> data4(new double[blocksize]);
          sort_indices<2,3,0,1,0,1,1,1>(data2.get()          , data3.get(), i0.size(), i1.size(), i2.size(), i3.size());
          sort_indices<0,3,2,1,0,1,1,1>(data2.get()+blocksize, data4.get(), i0.size(), i1.size(), i2.size(), i3.size());
          res->add_block(data3, i2, i3, i0, i1);
          res->add_block(data4, i0, i3, i2, i1);
        }
      }
    }
    }
  }

  // a r s t case
  {
    for (auto& i3 : active_) {
    for (auto& i2 : active_) {
    for (auto& i0 : active_) {
      assert(denom_->shalf_xxh());
      const size_t interm_size = denom_->shalf_xxh()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I2, const Index& I3) {
        unique_ptr<double[]> out(new double[I0.size()*I2.size()*I3.size()*interm_size]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
          for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
            for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
              copy_n(denom_->shalf_xxh()->element_ptr(0,j0-nclosed+nact*(j2-nclosed+nact*(j3-nclosed)) + i*nact*nact*nact),
                     interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i0, i2, i3);

      for (auto& i1 : virt_) {
        if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
        const size_t blocksize = norm->get_size(i2, i3, i0, i1);
        // data0 is the source area
        unique_ptr<double[]> data0 = norm->get_block(i2, i3, i0, i1);
        unique_ptr<double[]> data1(new double[blocksize]);
        // sort. Active indices run slower
        sort_indices<3,2,0,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
        // intermediate area
        unique_ptr<double[]> interm(new double[i1.size()*interm_size]);

        // move to orthogonal basis
        btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i1.size(), interm_size, i0.size()*i2.size()*i3.size(),
                                    1.0, data1.get(), i1.size(), transp.get(), interm_size, 0.0, interm.get(), i1.size());

        size_t iall = 0;
        for (int j123 = 0; j123 != interm_size; ++j123)
          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++iall)
            interm[iall] *= shift2 / (denom_->denom_xxh(j123) + eig_[j1] - e0_);

        // move back to non-orthogonal basis
        btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i1.size(), i0.size()*i2.size()*i3.size(), interm_size,
                                    1.0, interm.get(), i1.size(), transp.get(), interm_size, 0.0, data1.get(), i1.size());

        // sort back to the original order
        unique_ptr<double[]> data2(new double[blocksize]);
        sort_indices<2,3,1,0,0,1,1,1>(data1, data2, i1.size(), i0.size(), i2.size(), i3.size());
        res->add_block(data2, i2, i3, i0, i1);
      }
    }
    }
    }
  }

  // r i s t case
  {
    for (auto& i3 : active_) {
    for (auto& i1 : active_) {
    for (auto& i0 : active_) {
      assert(denom_->shalf_xhh());
      const size_t interm_size = denom_->shalf_xhh()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I1, const Index& I3) {
        unique_ptr<double[]> out(new double[I0.size()*I1.size()*I3.size()*interm_size]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
          for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
            for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
              copy_n(denom_->shalf_xhh()->element_ptr(0,j0-nclosed+nact*(j1-nclosed+nact*(j3-nclosed)) + i*nact*nact*nact),
                     interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<double[]> transp = create_transp(istate, i0, i1, i3);

      for (auto& i2 : closed_) {
        if (!norm->is_local(i2, i3, i0, i1) || !res->get_size(i2, i3, i0, i1)) continue;
        // if this block is not included in the current wave function, skip it
        const size_t blocksize = norm->get_size(i2, i3, i0, i1);
        // data0 is the source area
        unique_ptr<double[]> data0 = norm->get_block(i2, i3, i0, i1);
        unique_ptr<double[]> data1(new double[blocksize]);
        // sort. Active indices run slower
        sort_indices<0,2,3,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
        // intermediate area
        unique_ptr<double[]> interm(new double[i2.size()*interm_size]);

        // move to orthogonal basis
        btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i2.size(), interm_size, i0.size()*i1.size()*i3.size(),
                                    1.0, data1.get(), i2.size(), transp.get(), interm_size, 0.0, interm.get(), i2.size());

        size_t iall = 0;
        for (int j013 = 0; j013 != interm_size; ++j013)
          for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++iall)
            interm[iall] *= shift2 / (denom_->denom_xhh(j013) - eig_[j2] - e0_);

        // move back to non-orthogonal basis
        btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i2.size(), i0.size()*i1.size()*i3.size(), interm_size,
                                    1.0, interm.get(), i2.size(), transp.get(), interm_size, 0.0, data1.get(), i2.size());

        // sort back to the original order
        unique_ptr<double[]> data2(new double[blocksize]);
        sort_indices<0,3,1,2,0,1,1,1>(data1, data2, i2.size(), i0.size(), i1.size(), i3.size());
        res->add_block(data2, i2, i3, i0, i1);
      }
    }
    }
    }
  }
#endif

}

#endif
