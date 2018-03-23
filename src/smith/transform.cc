//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: transform.cc
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

#include <numeric>
#include <src/smith/moint.h>
#include <src/smith/spinfreebase.h>
#include <src/smith/smith_util.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template<typename DataType>
shared_ptr<Vector_<DataType>> SpinFreeMethod<DataType>::transform_to_orthogonal(shared_ptr<const MultiTensor_<DataType>> tensor) const {
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;

  const int nst = tensor->nref();
  shared_ptr<Vector_<DataType>> out;
  if (info_->sssr())
    out = make_shared<Vector_<DataType>>(size_all);
  else
    out = make_shared<Vector_<DataType>>(nst*size_all);

  int ioffset = 0;
  for (int ist = 0; ist != nst; ++ist) {
    if (tensor->at(ist)) {
      shared_ptr<Vector_<DataType>> vist = transform_to_orthogonal(tensor->at(ist), ist);
      copy_n(vist->data(), vist->size(), out->data()+ioffset);
      ioffset += vist->size();
    }
  }

  return out;
}



template<typename DataType>
shared_ptr<Vector_<DataType>> SpinFreeMethod<DataType>::transform_to_orthogonal(shared_ptr<const Tensor_<DataType>> tensor, const int istate) const {
  // number of orthogonal basis functions: use numbers of denom.
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;

  auto out = make_shared<Vector_<DataType>>(size_all);

  // a i b j case
  size_t ioffset = 0;
  {
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            unique_ptr<DataType[]> data0 = tensor->get_block(i0, i1, i2, i3);
            // not sure whether this is nessessary or not...
            if (is_same<DataType,double>::value) {
              const unique_ptr<DataType[]> data1 = tensor->get_block(i0, i3, i2, i1);
              sort_indices<0,3,2,1,2,12,1,12>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
            } else {
              blas::scale_n(0.25, data0.get(), tensor->get_size(i0, i1, i2, i3));
            }
            size_t iall = 0;
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                  for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                    const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                    (*out)[ioffset + jall] = data0[iall];
                  }
            }
          }
  }

  // a r b s case
  {
    ioffset += size_aibj;
    for (auto& i2 : active_)
      for (auto& i0 : active_) {
        const size_t interm_size = denom_->shalf_xx()->ndim();
        auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I2) {
          unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*interm_size]);
          for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
            for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
              copy_n(denom_->shalf_xx()->element_ptr(0,(j0-nclosed)+(j2-nclosed)*nact + i*nact*nact),
                     interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(istate, i0, i2);
        for (auto& i3 : virt_)
          for (auto& i1 : virt_) {
            const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> data0 = tensor->get_block(i0, i1, i2, i3);
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
            unique_ptr<DataType[]> interm(new DataType[i1.size()*i3.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i3.size(), i0.size()*i2.size(),
                                        1.0, transp.get(), interm_size, data1.get(), i0.size()*i2.size(), 0.0, interm.get(), interm_size);
            size_t iall = 0;
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                for (int j02 = 0; j02 != interm_size; ++j02, ++iall) {
                  const size_t jall = j02 + interm_size * (j1 + nvirt * j3);
                  (*out)[ioffset + jall] = interm[iall];
                }
              }
            }
          }
      }
  }

  // a r b i case
  {
    ioffset += size_arbs;
    for (auto& i0 : active_) {
      const size_t interm_size = denom_->shalf_x()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0) {
        unique_ptr<DataType[]> out(new DataType[I0.size()*interm_size]);
        for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
          copy_n(denom_->shalf_x()->element_ptr(0,j0-nclosed + i*nact), interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<DataType[]> transp = create_transp(istate, i0);
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_) {
            const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
            unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data2(new DataType[blocksize]);
            sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
            if (is_same<DataType,double>::value) {
              const unique_ptr<DataType[]> data1 = tensor->get_block(i2, i1, i0, i3);
              sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());
            } else {
              blas::scale_n(0.5, data2.get(), blocksize);
            }
            // move to orthogonal basis
            unique_ptr<DataType[]> interm(new DataType[i1.size()*i2.size()*i3.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i2.size()*i3.size(), i0.size(),
                                        1.0, transp.get(), interm_size, data2.get(), i0.size(), 0.0, interm.get(), interm_size);

            size_t iall = 0;
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                  for (int j0 = 0; j0 != interm_size; ++j0, ++iall) {
                    const size_t jall = j0 + interm_size * (j1 + nvirt * (j2 + nclo * j3));
                    (*out)[ioffset + jall] = interm[iall];
                  }
          }
    }
  }

  // a i r j case
  {
    ioffset += size_arbi;
    for (auto& i3 : active_) {
      const size_t interm_size = denom_->shalf_h()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I3) {
        unique_ptr<DataType[]> out(new DataType[I3.size()*interm_size]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
          copy_n(denom_->shalf_h()->element_ptr(0,j3-nclosed + i*nact), interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<DataType[]> transp = create_transp(istate, i3);
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
            unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data2(new DataType[blocksize]);
            sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
            if (is_same<DataType,double>::value) {
              const unique_ptr<DataType[]> data1 = tensor->get_block(i0, i3, i2, i1);
              sort_indices<0,3,2,1,2,3,1,3>(data1, data2, i0.size(), i3.size(), i2.size(), i1.size());
            } else {
              blas::scale_n(0.5, data2.get(), blocksize);
            }
            unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*i2.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size()*i2.size(), interm_size, i3.size(),
                                        1.0, data2.get(), i0.size()*i1.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size()*i2.size());
            size_t iall = 0;
            for (int j3 = 0; j3 != interm_size; ++j3)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                  for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                    const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                    (*out)[ioffset + jall] = interm[iall];
                  }
          }
    }
  }

  // r i s j case
  {
    ioffset += size_airj;
    for (auto& i3 : active_)
      for (auto& i1 : active_) {
        const size_t interm_size = denom_->shalf_hh()->ndim();
        auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I1, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I1.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
              copy_n(denom_->shalf_hh()->element_ptr(0,(j1-nclosed)+(j3-nclosed)*nact + i*nact*nact),
                     interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(istate, i1, i3);
        for (auto& i2 : closed_)
          for (auto& i0 : closed_) {
            const size_t blocksize = tensor->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> data0 = tensor->get_block(i0, i1, i2, i3);
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
            unique_ptr<DataType[]> interm(new DataType[i0.size()*i2.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i2.size(), interm_size, i1.size()*i3.size(),
                                        1.0, data1.get(), i0.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i2.size());
            size_t iall = 0;
            for (int j13 = 0; j13 != interm_size; ++j13)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                  const size_t jall = j0 + nclo * (j2 + nclo * j13);
                  (*out)[ioffset + jall] = interm[iall];
                }
          }
      }
  }

  // a i r s & a r s i case
  // TODO implement complex case
  {
    ioffset += size_risj;
    for (auto& i3 : active_)
      for (auto& i2 : active_) {
        const size_t interm_size = denom_->shalf_xh()->ndim();
        auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*interm_size*2]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclosed)+(j3-nclosed)*nact + 2*i*nact*nact),
                     interm_size, out.get()+interm_size*k);
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclosed)+(j3-nclosed)*nact + (2*i+1)*nact*nact),
                     interm_size, out.get()+interm_size*(k+I2.size()*I3.size()));
            }
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(istate, i2, i3);
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
            unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data1 = tensor->get_block(i0, i3, i2, i1);
            unique_ptr<DataType[]> data2(new DataType[blocksize*2]);
            sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
            sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
            unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size(), interm_size, i2.size()*i3.size()*2,
                                        1.0, data2.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size());
            size_t iall = 0;
            for (int j23 = 0; j23 != interm_size; ++j23)
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                  const size_t jall = j0 + nclo * (j1 + nvirt * j23);
                  (*out)[ioffset + jall] = interm[iall];
                }
          }
      }
  }

  // a r s t case
  {
    ioffset += size_airs;
    for (auto& i3 : active_)
      for (auto& i2 : active_)
        for (auto& i0 : active_) {
          const size_t interm_size = denom_->shalf_xxh()->ndim();
          auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I2, const Index& I3) {
            unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*I3.size()*interm_size]);
            for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
              for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
                for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                  copy_n(denom_->shalf_xxh()->element_ptr(0,j0-nclosed+nact*(j2-nclosed+nact*(j3-nclosed)) + i*nact*nact*nact),
                         interm_size, out.get()+interm_size*k);
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(istate, i0, i2, i3);
          for (auto& i1 : virt_) {
            const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
            unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<3,2,0,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
            unique_ptr<DataType[]> interm(new DataType[i1.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i1.size(), interm_size, i0.size()*i2.size()*i3.size(),
                                        1.0, data1.get(), i1.size(), transp.get(), interm_size, 0.0, interm.get(), i1.size());
            size_t iall = 0;
            for (int j023 = 0; j023 != interm_size; ++j023)
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1, ++iall) {
                const size_t jall = j1 + nvirt * j023;
                (*out)[ioffset + jall] = interm[iall];
              }
          }
        }
  }

  // r i s t case
  {
    ioffset += size_arst;
    for (auto& i3 : active_)
      for (auto& i1 : active_)
        for (auto& i0 : active_) {
          const size_t interm_size = denom_->shalf_xhh()->ndim();
          auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I1, const Index& I3) {
            unique_ptr<DataType[]> out(new DataType[I0.size()*I1.size()*I3.size()*interm_size]);
            for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
              for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
                for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                  copy_n(denom_->shalf_xhh()->element_ptr(0,j0-nclosed+nact*(j1-nclosed+nact*(j3-nclosed)) + i*nact*nact*nact),
                         interm_size, out.get()+interm_size*k);
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(istate, i0, i1, i3);
          for (auto& i2 : closed_) {
            const size_t blocksize = tensor->get_size(i2, i3, i0, i1);
            unique_ptr<DataType[]> data0 = tensor->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<0,2,3,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
            unique_ptr<DataType[]> interm(new DataType[i2.size()*interm_size]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i2.size(), interm_size, i0.size()*i1.size()*i3.size(),
                                        1.0, data1.get(), i2.size(), transp.get(), interm_size, 0.0, interm.get(), i2.size());
            size_t iall = 0;
            for (int j013 = 0; j013 != interm_size; ++j013)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2, ++iall) {
                const size_t jall = j2 + nclo*j013;
                (*out)[ioffset + jall] = interm[iall];
              }
          }
        }
    ioffset += size_rist;
  }

  return out;
}


template<typename DataType>
shared_ptr<MultiTensor_<DataType>> SpinFreeMethod<DataType>::transform_to_redundant_residual(shared_ptr<const Vector_<DataType>> vector, const int nstates, const int istate) const {
  // not urgent then transform_to_redundant
  auto out = make_shared<MultiTensor_<DataType>>(nstates);

  return out;
}


template<typename DataType>
shared_ptr<Tensor_<DataType>> SpinFreeMethod<DataType>::transform_to_redundant_residual_v(shared_ptr<const Vector_<DataType>> vector, const int istate) const {
  auto out = init_residual();

  return out;
}


template<typename DataType>
shared_ptr<MultiTensor_<DataType>> SpinFreeMethod<DataType>::transform_to_redundant_amplitude(shared_ptr<const Vector_<DataType>> vector, const int nstates, const int istate) const {
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;
  auto out = make_shared<MultiTensor_<DataType>>(nstates);

  size_t ioffset = 0;
  for (int i = 0; i != nstates; ++i) {
    if (!info_->sssr() || i == istate) {
      auto vi = make_shared<Vector_<DataType>>(size_all);
      copy_n(vector->data()+ioffset, size_all, vi->data());
      (*out)[i] = transform_to_redundant_amplitude_v(vi, i);
      ioffset += size_all;
    }
  }

  return out;
}


template<typename DataType>
shared_ptr<Tensor_<DataType>> SpinFreeMethod<DataType>::transform_to_redundant_amplitude_v(shared_ptr<const Vector_<DataType>> vector, const int istate) const {
  // to be implemented
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

  auto out = init_amplitude();

  size_t ioffset = 0;

  // a i b j case
  {
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            size_t iall = 0;
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                  for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                    const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));

                    data0[iall] = (*vector)[ioffset + jall];
                  }
            }
            out->put_block(data0, i0, i1, i2, i3);
          }
  }

  // a r b s case
  {
    ioffset += size_aibj;
    for (auto& i2 : active_)
      for (auto& i0 : active_) {
        const size_t interm_size = denom_->shalf_xx()->ndim();
        auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I2) {
          unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*interm_size]);
          for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
            for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
              copy_n(denom_->shalf_xx()->element_ptr(0,(j0-nclosed)+(j2-nclosed)*nact + i*nact*nact),
                     interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(istate, i0, i2);
        for (auto& i3 : virt_)
          for (auto& i1 : virt_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> interm(new DataType[i1.size()*i3.size()*interm_size]);

            size_t iall = 0;
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                for (int j02 = 0; j02 != interm_size; ++j02, ++iall) {
                  const size_t jall = j02 + interm_size * (j1 + nvirt * j3);
                  interm[iall] = (*vector)[ioffset + jall];
                }
              }
            }

            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), interm_size,
                                        0.5, transp.get(), interm_size, interm.get(), interm_size, 0.0, data0.get(), i0.size()*i2.size());
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
            out->put_block(data1, i0, i1, i2, i3);
          }
      }
  }

  // a r b i case
  {
    ioffset += size_arbs;
    for (auto& i0 : active_) {
      const size_t interm_size = denom_->shalf_x()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0) {
        unique_ptr<DataType[]> out(new DataType[I0.size()*interm_size]);
        for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
          copy_n(denom_->shalf_x()->element_ptr(0,j0-nclosed + i*nact), interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<DataType[]> transp = create_transp(istate, i0);
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> interm(new DataType[i1.size()*i2.size()*i3.size()*interm_size]);

            size_t iall = 0;
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                  for (int j0 = 0; j0 != interm_size; ++j0, ++iall) {
                    const size_t jall = j0 + interm_size * (j1 + nvirt * (j2 + nclo * j3));
                    interm[iall] = (*vector)[ioffset + jall];
                  }

            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size(), i1.size()*i2.size()*i3.size(), interm_size,
                                        1.0, transp.get(), interm_size, interm.get(), interm_size, 0.0, data0.get(), i0.size());

            out->put_block(data0, i0, i1, i2, i3);
          }
    }
  }

  // a i r j case
  {
    ioffset += size_arbi;
    for (auto& i3 : active_) {
      const size_t interm_size = denom_->shalf_h()->ndim();
      auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I3) {
        unique_ptr<DataType[]> out(new DataType[I3.size()*interm_size]);
        for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
          copy_n(denom_->shalf_h()->element_ptr(0,j3-nclosed + i*nact), interm_size, out.get()+interm_size*k);
        return move(out);
      };
      unique_ptr<DataType[]> transp = create_transp(istate, i3);
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);

            unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*i2.size()*interm_size]);
            size_t iall = 0;
            for (int j3 = 0; j3 != interm_size; ++j3)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                  for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                    const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                    interm[iall] = (*vector)[ioffset + jall];
                  }

            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size()*i2.size(), i3.size(), interm_size,
                                        1.0, interm.get(), i0.size()*i1.size()*i2.size(), transp.get(), interm_size, 0.0, data0.get(), i0.size()*i1.size()*i2.size());

            out->put_block(data0, i0, i1, i2, i3);
          }
    }
  }

  // r i s j case
  {
    ioffset += size_airj;
    for (auto& i3 : active_)
      for (auto& i1 : active_) {
        const size_t interm_size = denom_->shalf_hh()->ndim();
        auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I1, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I1.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
              copy_n(denom_->shalf_hh()->element_ptr(0,(j1-nclosed)+(j3-nclosed)*nact + i*nact*nact),
                     interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(istate, i1, i3);
        for (auto& i2 : closed_)
          for (auto& i0 : closed_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> interm(new DataType[i0.size()*i2.size()*interm_size]);
            size_t iall = 0;
            for (int j13 = 0; j13 != interm_size; ++j13)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                  const size_t jall = j0 + nclo * (j2 + nclo * j13);
                  interm[iall] = (*vector)[ioffset + jall];
                }
            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), interm_size,
                                        0.5, interm.get(), i0.size()*i2.size(), transp.get(), interm_size, 0.0, data0.get(), i0.size()*i2.size());
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
            out->put_block(data1, i0, i1, i2, i3);
          }
      }
  }

  // a i r s & a r s i case
  // TODO implement complex case
  {
    ioffset += size_risj;
    for (auto& i3 : active_)
      for (auto& i2 : active_) {
        const size_t interm_size = denom_->shalf_xh()->ndim();
        auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*interm_size*2]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclosed)+(j3-nclosed)*nact + 2*i*nact*nact),
                     interm_size, out.get()+interm_size*k);
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclosed)+(j3-nclosed)*nact + (2*i+1)*nact*nact),
                     interm_size, out.get()+interm_size*(k+I2.size()*I3.size()));
            }
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(istate, i2, i3);
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*interm_size]);
            size_t iall = 0;
            for (int j23 = 0; j23 != interm_size; ++j23)
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0, ++iall) {
                  const size_t jall = j0 + nclo * (j1 + nvirt * j23);
                  interm[iall] = (*vector)[ioffset + jall];
               }
            unique_ptr<DataType[]> data0(new DataType[blocksize*2]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2.size()*i3.size()*2, interm_size,
                                        1.0, interm.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, data0.get(), i0.size()*i1.size());
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            unique_ptr<DataType[]> data2(new DataType[blocksize]);
            copy_n(data0.get(), blocksize, data1.get());
            sort_indices<2,1,0,3,0,1,1,1>(data0.get()+blocksize, data2.get(), i0.size(), i1.size(), i2.size(), i3.size());
            out->put_block(data1, i0, i1, i2, i3);
            out->put_block(data2, i2, i1, i0, i3);
          }
      }
  }

  // a r s t case
  {
    ioffset += size_airs;
    for (auto& i3 : active_)
      for (auto& i2 : active_)
        for (auto& i0 : active_) {
          const size_t interm_size = denom_->shalf_xxh()->ndim();
          auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I2, const Index& I3) {
            unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*I3.size()*interm_size]);
            for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
              for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
                for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                  copy_n(denom_->shalf_xxh()->element_ptr(0,j0-nclosed+nact*(j2-nclosed+nact*(j3-nclosed)) + i*nact*nact*nact),
                         interm_size, out.get()+interm_size*k);
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(istate, i0, i2, i3);
          for (auto& i1 : virt_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> interm(new DataType[i1.size()*interm_size]);
            size_t iall = 0;
            for (int j023 = 0; j023 != interm_size; ++j023)
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1, ++iall) {
                const size_t jall = j1 + nvirt * j023;
                interm[iall] = (*vector)[ioffset + jall];
              }
            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i1.size(), i0.size()*i2.size()*i3.size(), interm_size,
                                        1.0, interm.get(), i1.size(), transp.get(), interm_size, 0.0, data0.get(), i1.size());
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<1,0,2,3,0,1,1,1>(data0, data1, i1.size(), i0.size(), i2.size(), i3.size());
            out->put_block(data1, i0, i1, i2, i3);
          }
        }
  }

  // r i s t case
  {
    ioffset += size_arst;
    for (auto& i3 : active_)
      for (auto& i1 : active_)
        for (auto& i0 : active_) {
          const size_t interm_size = denom_->shalf_xhh()->ndim();
          auto create_transp = [&nclosed,&nact,&interm_size, this](const int i, const Index& I0, const Index& I1, const Index& I3) {
            unique_ptr<DataType[]> out(new DataType[I0.size()*I1.size()*I3.size()*interm_size]);
            for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
              for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
                for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                  copy_n(denom_->shalf_xhh()->element_ptr(0,j0-nclosed+nact*(j1-nclosed+nact*(j3-nclosed)) + i*nact*nact*nact),
                         interm_size, out.get()+interm_size*k);
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(istate, i0, i1, i3);
          for (auto& i2 : closed_) {
            const size_t blocksize = out->get_size(i0, i1, i2, i3);
            unique_ptr<DataType[]> interm(new DataType[i2.size()*interm_size]);
            size_t iall = 0;
            for (int j013 = 0; j013 != interm_size; ++j013)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2, ++iall) {
                const size_t jall = j2 + nclo*j013;
                interm[iall] = (*vector)[ioffset + jall];
              }
            unique_ptr<DataType[]> data0(new DataType[blocksize]);
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i2.size(), i0.size()*i1.size()*i3.size(), interm_size,
                                        1.0, interm.get(), i2.size(), transp.get(), interm_size, 0.0, data0.get(), i2.size());
            unique_ptr<DataType[]> data1(new DataType[blocksize]);
            sort_indices<1,2,0,3,0,1,1,1>(data0, data1, i2.size(), i0.size(), i1.size(), i3.size());
            out->put_block(data1, i0, i1, i2, i3);
          }
        }

    ioffset += size_rist;
  }

  return out;
}


template<typename DataType>
void SpinFreeMethod<DataType>::update_amplitude_orthogonal(shared_ptr<Vector_<DataType>> amplitude, shared_ptr<const Vector_<DataType>> residual, const int nstates, const int istate) const {
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;
  const double shift2 = info_->shift() * info_->shift();

  size_t ioffset = 0;

  for (int ist = 0; ist != nstates; ++ist) {
    if (!info_->sssr() || ist == istate) {

      // a i b j case
      {
        const double e0loc = e0all_[ist] - e0_;
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0 : closed_) {
                for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                      for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                        const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                        const double denom = e0loc - eig_[j0+ncore] - eig_[j2+ncore] + eig_[j3+nocc] + eig_[j1+nocc];
                        if (info_->shift_imag()) {
                          (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                        } else {
                          (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                        }
                      }
                }
              }
      }

      // a r b s case
      {
        ioffset += size_aibj;
        for (auto& i2 : active_)
          for (auto& i0 : active_) {
            const size_t interm_size = denom_->shalf_xx()->ndim();
            for (auto& i3 : virt_)
              for (auto& i1 : virt_) {
                for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                    for (int j02 = 0; j02 != interm_size; ++j02) {
                      const size_t jall = j02 + interm_size * (j1 + nvirt * j3);
                      const double denom = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j02) - e0_;
                      if (info_->shift_imag()) {
                        (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                      } else {
                        (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                      }
                    }
                  }
                }
              }
          }
      }

      // a r b i case
      {
        ioffset += size_arbs;
        for (auto& i0 : active_) {
          const size_t interm_size = denom_->shalf_x()->ndim();
          for (auto& i3 : virt_)
            for (auto& i2 : closed_)
              for (auto& i1 : virt_) {
                for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                      for (int j0 = 0; j0 != interm_size; ++j0) {
                        const size_t jall = j0 + interm_size * (j1 + nvirt * (j2 + nclo * j3));
                        const double denom = eig_[j1+nocc] + eig_[j3+nocc] - eig_[j2+ncore] + denom_->denom_x(j0) - e0_;
                        if (info_->shift_imag()) {
                          (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                        } else {
                          (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                        }
                      }
              }
        }
      }

      // a i r j case
      {
        ioffset += size_arbi;
        for (auto& i3 : active_) {
          const size_t interm_size = denom_->shalf_h()->ndim();
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0 : closed_) {
                for (int j3 = 0; j3 != interm_size; ++j3)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                      for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                        const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                        const double denom = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j3) - e0_;
                        if (info_->shift_imag()) {
                          (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                        } else {
                          (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                        }
                      }
              }
        }
      }

      // r i s j case
      {
        ioffset += size_airj;
        for (auto& i3 : active_)
          for (auto& i1 : active_) {
            const size_t interm_size = denom_->shalf_hh()->ndim();
            for (auto& i2 : closed_)
              for (auto& i0 : closed_) {
                for (int j13 = 0; j13 != interm_size; ++j13)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                      const size_t jall = j0 + nclo * (j2 + nclo * j13);
                      const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_hh(j13) - e0_;
                      if (info_->shift_imag()) {
                        (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                      } else {
                        (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                      }
                    }
              }
          }
      }

      // a i r s & a r s i case
      // TODO implement complex case
      {
        ioffset += size_risj;
        for (auto& i3 : active_)
          for (auto& i2 : active_) {
            const size_t interm_size = denom_->shalf_xh()->ndim();
            for (auto& i1 : virt_)
              for (auto& i0 : closed_) {
                for (int j23 = 0; j23 != interm_size; ++j23)
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                    for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                      const size_t jall = j0 + nclo * (j1 + nvirt * j23);
                      const double denom = eig_[j1+nocc] - eig_[j0+ncore] + denom_->denom_xh(j23) - e0_;
                      if (info_->shift_imag()) {
                        (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                      } else {
                        (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                      }
                   }
              }
          }
      }

      // a r s t case
      {
        ioffset += size_airs;
        for (auto& i3 : active_)
          for (auto& i2 : active_)
            for (auto& i0 : active_) {
              const size_t interm_size = denom_->shalf_xxh()->ndim();
              for (auto& i1 : virt_) {
                for (int j023 = 0; j023 != interm_size; ++j023)
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                    const size_t jall = j1 + nvirt * j023;
                    const double denom = eig_[j1+nocc] + denom_->denom_xxh(j023) - e0_;
                    if (info_->shift_imag()) {
                      (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                    } else {
                      (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                    }
                  }
              }
            }
      }

      // r i s t case
      {
        ioffset += size_arst;
        for (auto& i3 : active_)
          for (auto& i1 : active_)
            for (auto& i0 : active_) {
              const size_t interm_size = denom_->shalf_xhh()->ndim();
              for (auto& i2 : closed_) {
                for (int j013 = 0; j013 != interm_size; ++j013)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2) {
                    const size_t jall = j2 + nclo*j013;
                    const double denom = - eig_[j2+ncore] + denom_->denom_xhh(j013) - e0_;
                    if (info_->shift_imag()) {
                      (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] * denom / (denom * denom + shift2);
                    } else {
                      (*amplitude)[ioffset + jall] -= (*residual)[ioffset + jall] / denom;
                    }
                  }
              }
            }
        ioffset += size_rist;
      }

    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explicit instantiation at the end of the file
template class SpinFreeMethod<double>;
template class SpinFreeMethod<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
