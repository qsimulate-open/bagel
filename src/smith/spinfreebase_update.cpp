//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase_denom.cpp
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#ifdef SPINFREEMETHOD_DETAIL

template<typename DataType>
void SpinFreeMethod<DataType>::update_amplitude(shared_ptr<Tensor_<DataType>> t, shared_ptr<const Tensor_<DataType>> r) const {
  shared_ptr<MultiTensor_<DataType>> tt
    = make_shared<MultiTensor_<DataType>>(vector<DataType>{0.0}, vector<shared_ptr<Tensor_<DataType>>>{t});

  auto r0 = const_pointer_cast<Tensor_<DataType>>(r);
  shared_ptr<const MultiTensor_<DataType>> rr
    = make_shared<MultiTensor_<DataType>>(vector<DataType>{0.0}, vector<shared_ptr<Tensor_<DataType>>>{r0});

  update_amplitude(tt, rr);
}


template<typename DataType>
void SpinFreeMethod<DataType>::update_amplitude(shared_ptr<MultiTensor_<DataType>> t, shared_ptr<const MultiTensor_<DataType>> r) const {

  if (t->nref() != r->nref())
    throw logic_error("something is wrong. SpinFreeMethod::update_amplitude");

  const int nst = t->nref();
  const int fac2 = is_same<DataType,double>::value ? 1.0 : 2.0;

  for (int ist = 0; ist != nst; ++ist) {
    t->fac(ist) = 0.0;
    for (auto& i3 : virt_) {
      for (auto& i2 : closed_) {
        for (auto& i1 : virt_) {
          for (auto& i0 : closed_) {
            // if this block is not included in the current wave function, skip it
            if (!r->at(ist)->get_size_alloc(i0, i1, i2, i3)) continue;
            unique_ptr<DataType[]>       data0 = r->at(ist)->get_block(i0, i1, i2, i3);

            // this is an inverse of the overlap.
            if (is_same<DataType,double>::value) {
              const unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i0, i3, i2, i1);
              sort_indices<0,3,2,1,2,12,1,12>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
            } else {
              blas::scale_n(0.25, data0.get(), r->at(ist)->get_size_alloc(i0, i1, i2, i3));
            }
            size_t iall = 0;
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    // note that e0 is cancelled by another term
                    data0[iall] /= eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1];
            t->at(ist)->add_block(data0, i0, i1, i2, i3);
          }
        }
      }
    }

    // not the best structure, but I am assuming that this does not take too much time...
    for (int jst = 0; jst != nst; ++jst) {

      for (auto& i2 : active_) {
        for (auto& i0 : active_) {
          // trans is the transformation matrix
          assert(denom_->shalf_xx());
          const size_t interm_size = denom_->shalf_xx()->ndim();
          const int nact = info_->nact() * fac2;
          const int nclo = info_->nclosed() * fac2;
          auto create_transp = [&,this](const int i) {
            unique_ptr<DataType[]> out(new DataType[i0.size()*i2.size()*interm_size]);
            for (int j2 = i2.offset(), k = 0; j2 != i2.offset()+i2.size(); ++j2)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                copy_n(denom_->shalf_xx()->element_ptr(0,(j0-nclo)+(j2-nclo)*nact + i*nact*nact),
                       interm_size, out.get()+interm_size*k);
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(ist);
          unique_ptr<DataType[]> transp2 = create_transp(jst);

          for (auto& i3 : virt_) {
            for (auto& i1 : virt_) {
              // if this block is not included in the current wave function, skip it
              if (!r->at(ist)->get_size_alloc(i0, i1, i2, i3)) continue;
              // data0 is the source area
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i0, i1, i2, i3);
              unique_ptr<DataType[]> data1(new DataType[r->at(ist)->get_size(i0, i1, i2, i3)]);
              // sort. Active indices run faster
              if (is_same<DataType,double>::value)
                sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              else
                sort_indices<0,2,1,3,0,1,2,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i1.size()*i3.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i3.size(), i0.size()*i2.size(),
                                          1.0, transp.get(), interm_size, data1.get(), i0.size()*i2.size(), 0.0, interm.get(), interm_size);

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j02 = 0; j02 != interm_size; ++j02, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom_xx(j02) + eig_[j3] + eig_[j1]));

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              // TODO check for complex cases
              btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), interm_size,
                                          0.5, transp2.get(), interm_size, interm.get(), interm_size, 0.0, data0.get(), i0.size()*i2.size());

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
              t->at(jst)->add_block(data1, i0, i1, i2, i3);
            }
          }
        }
      }

      for (auto& i0 : active_) {
        // trans is the transformation matrix
        assert(denom_->shalf_x());
        const size_t interm_size = denom_->shalf_x()->ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&,this](const int i) {
          unique_ptr<DataType[]> out(new DataType[i0.size()*interm_size]);
          for (int j0 = i0.offset(), k = 0; j0 != i0.offset()+i0.size(); ++j0, ++k)
            copy_n(denom_->shalf_x()->element_ptr(0,j0-nclo + i*nact), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ist);
        unique_ptr<DataType[]> transp2 = create_transp(jst);

        for (auto& i3 : virt_) {
          for (auto& i2 : closed_) {
            for (auto& i1 : virt_) {
              if (!r->at(ist)->get_size_alloc(i2, i3, i0, i1)) continue;
              unique_ptr<DataType[]>       data0 = r->at(ist)->get_block(i2, i3, i0, i1);
              unique_ptr<DataType[]> data2(new DataType[r->at(ist)->get_size(i2, i3, i0, i1)]);
              sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
              if (is_same<DataType,double>::value) {
                assert(r->at(ist)->get_size_alloc(i2, i1, i0, i3));
                const unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i2, i1, i0, i3);
                sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());
              } else {
                blas::scale_n(0.5, data2.get(), r->at(ist)->get_size(i2, i3, i0, i1));
              }

              // move to orthogonal basis
              unique_ptr<DataType[]> interm(new DataType[i1.size()*i2.size()*i3.size()*interm_size]);
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i2.size()*i3.size(), i0.size(),
                                          1.0, transp.get(), interm_size, data2.get(), i0.size(), 0.0, interm.get(), interm_size);

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = 0; j0 != interm_size; ++j0, ++iall)
                      interm[iall] /= min(-0.1, e0_ - (denom_->denom_x(j0) + eig_[j3] - eig_[j2] + eig_[j1]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0.size(), i1.size()*i2.size()*i3.size(), interm_size,
                                          1.0, transp2.get(), interm_size, interm.get(), interm_size, 0.0, data2.get(), i0.size());

              t->at(jst)->add_block(data2, i0, i1, i2, i3);
            }
          }
        }
      }

      for (auto& i3 : active_) {
        // trans is the transformation matrix
        assert(denom_->shalf_h());
        const size_t interm_size = denom_->shalf_x()->ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&,this](const int i) {
          unique_ptr<DataType[]> out(new DataType[i3.size()*interm_size]);
          for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3, ++k)
            copy_n(denom_->shalf_h()->element_ptr(0,j3-nclo + i*nact), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ist);
        unique_ptr<DataType[]> transp2 = create_transp(jst);

        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              if (!r->at(ist)->get_size_alloc(i2, i3, i0, i1)) continue;
              assert(r->at(ist)->get_size_alloc(i0, i3, i2, i1));
              unique_ptr<DataType[]>       data0 = r->at(ist)->get_block(i2, i3, i0, i1);
              const unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i0, i3, i2, i1);
              unique_ptr<DataType[]> data2(new DataType[r->at(ist)->get_size(i2, i3, i0, i1)]);
              sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
              sort_indices<0,3,2,1,2,3,1,3>(data1, data2, i0.size(), i3.size(), i2.size(), i1.size());
              unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*i2.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasConjTrans, i0.size()*i1.size()*i2.size(), interm_size, i3.size(),
                                          1.0, data2.get(), i0.size()*i1.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size()*i2.size());

              size_t iall = 0;
              for (int j3 = 0; j3 != interm_size; ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                      interm[iall] /= min(-0.1, e0_ - (denom_->denom_h(j3) - eig_[j2] + eig_[j1] - eig_[j0]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size()*i2.size(), i3.size(), interm_size,
                                          1.0, interm.get(), i0.size()*i1.size()*i2.size(), transp2.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size()*i2.size());

              t->at(jst)->add_block(data2, i0, i1, i2, i3);
            }
          }
        }
      }

      for (auto& i3 : active_) {
        for (auto& i1 : active_) {
          assert(denom_->shalf_hh());
          const size_t interm_size = denom_->shalf_hh()->ndim();
          const int nact = info_->nact() * fac2;
          const int nclo = info_->nclosed() * fac2;
          auto create_transp = [&,this](const int i) {
            unique_ptr<DataType[]> out(new DataType[i1.size()*i3.size()*interm_size]);
            for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++k)
                copy_n(denom_->shalf_hh()->element_ptr(0,(j1-nclo)+(j3-nclo)*nact + i*nact*nact),
                       interm_size, out.get()+interm_size*k);
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(ist);
          unique_ptr<DataType[]> transp2 = create_transp(jst);

          for (auto& i2 : closed_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              if (!r->at(ist)->get_size_alloc(i0, i1, i2, i3)) continue;
              // data0 is the source area
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i0, i1, i2, i3);
              unique_ptr<DataType[]> data1(new DataType[r->at(ist)->get_size(i0, i1, i2, i3)]);
              // sort. Active indices run slower
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i0.size()*i2.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasConjTrans, i0.size()*i2.size(), interm_size, i1.size()*i3.size(),
                                          1.0, data1.get(), i0.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i2.size());

              size_t iall = 0;
              for (int j13 = 0; j13 != interm_size; ++j13)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom_hh(j13) - eig_[j2] - eig_[j0]));

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i2.size(), i1.size()*i3.size(), interm_size,
                                          0.5, interm.get(), i0.size()*i2.size(), transp2.get(), interm_size, 0.0, data0.get(), i0.size()*i2.size());

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
              t->at(jst)->add_block(data1, i0, i1, i2, i3);
            }
          }
        }
      }

      for (auto& i3 : active_) {
        for (auto& i2 : active_) {
          assert(denom_->shalf_xh());
          const size_t interm_size = denom_->shalf_xh()->ndim();
          const int nact = info_->nact() * fac2;
          const int nclo = info_->nclosed() * fac2;
          auto create_transp = [&,this](const int i) {
            unique_ptr<DataType[]> out(new DataType[i2.size()*i3.size()*interm_size*2]);
            for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++k) {
                copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclo)+(j3-nclo)*nact + 2*i*nact*nact),
                       interm_size, out.get()+interm_size*k);
                copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclo)+(j3-nclo)*nact + (2*i+1)*nact*nact),
                       interm_size, out.get()+interm_size*(k+i2.size()*i3.size()));
              }
            return move(out);
          };
          unique_ptr<DataType[]> transp = create_transp(ist);
          unique_ptr<DataType[]> transp2 = create_transp(jst);

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size_alloc(i2, i3, i0, i1);
              if (!blocksize) continue;
              assert(blocksize == r->at(ist)->get_size_alloc(i0, i3, i2, i1));
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
              unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i0, i3, i2, i1);

              unique_ptr<DataType[]> data2(new DataType[blocksize*2]);
              // sort. Active indices run slower
              sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
              sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasConjTrans, i0.size()*i1.size(), interm_size, i2.size()*i3.size()*2,
                                          1.0, data2.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size());

              size_t iall = 0;
              for (int j23 = 0; j23 != interm_size; ++j23)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom_xh(j23) + eig_[j1] - eig_[j0]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2.size()*i3.size()*2, interm_size,
                                          1.0, interm.get(), i0.size()*i1.size(), transp2.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size());

              // sort back to the original order
              copy_n(data2.get(), blocksize, data0.get());
              sort_indices<2,1,0,3,0,1,1,1>(data2.get()+blocksize, data1.get(), i0.size(), i1.size(), i2.size(), i3.size());
              t->at(jst)->add_block(data0, i0, i1, i2, i3);
              t->at(jst)->add_block(data1, i2, i1, i0, i3);
            }
          }
        }
      }

      for (auto& i3 : active_) {
        for (auto& i2 : active_) {
          for (auto& i0 : active_) {
            assert(denom_->shalf_xhh());
            const size_t interm_size = denom_->shalf_xhh()->ndim();
            const int nact = info_->nact() * fac2;
            const int nclo = info_->nclosed() * fac2;
            auto create_transp = [&,this](const int i) {
              unique_ptr<DataType[]> out(new DataType[i0.size()*i2.size()*i3.size()*interm_size]);
              for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                    copy_n(denom_->shalf_xhh()->element_ptr(0,j0-nclo+nact*(j2-nclo+nact*(j3-nclo)) + i*nact*nact*nact),
                           interm_size, out.get()+interm_size*k);
              return move(out);
            };
            unique_ptr<DataType[]> transp = create_transp(ist);
            unique_ptr<DataType[]> transp2 = create_transp(jst);

            for (auto& i1 : virt_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size_alloc(i2, i3, i0, i1);
              if (!blocksize) continue;
              // data0 is the source area
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
              unique_ptr<DataType[]> data1(new DataType[blocksize]);
              // sort. Active indices run slower
              sort_indices<3,2,0,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i1.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasConjTrans, i1.size(), interm_size, i0.size()*i2.size()*i3.size(),
                                          1.0, data1.get(), i1.size(), transp.get(), interm_size, 0.0, interm.get(), i1.size());

              size_t iall = 0;
              for (int j123 = 0; j123 != interm_size; ++j123)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++iall)
                  interm[iall] /= min(-0.1, e0_ - (denom_->denom_xhh(j123) + eig_[j1]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i1.size(), i0.size()*i2.size()*i3.size(), interm_size,
                                          1.0, interm.get(), i1.size(), transp2.get(), interm_size, 0.0, data0.get(), i1.size());

              // sort back to the original order
              sort_indices<1,0,2,3,0,1,1,1>(data0, data1, i1.size(), i0.size(), i2.size(), i3.size());
              t->at(jst)->add_block(data1, i0, i1, i2, i3);
            }
          }
        }
      }

      for (auto& i3 : active_) {
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            assert(denom_->shalf_xxh());
            const size_t interm_size = denom_->shalf_xxh()->ndim();
            const int nact = info_->nact() * fac2;
            const int nclo = info_->nclosed() * fac2;
            auto create_transp = [&,this](const int i) {
              unique_ptr<DataType[]> out(new DataType[i0.size()*i1.size()*i3.size()*interm_size]);
              for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                    copy_n(denom_->shalf_xxh()->element_ptr(0,j0-nclo+nact*(j1-nclo+nact*(j3-nclo)) + i*nact*nact*nact),
                           interm_size, out.get()+interm_size*k);
              return move(out);
            };
            unique_ptr<DataType[]> transp = create_transp(ist);
            unique_ptr<DataType[]> transp2 = create_transp(jst);

            for (auto& i2 : closed_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size_alloc(i2, i3, i0, i1);
              if (!blocksize) continue;
              // data0 is the source area
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
              unique_ptr<DataType[]> data1(new DataType[blocksize]);
              // sort. Active indices run slower
              sort_indices<0,2,3,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i2.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasConjTrans, i2.size(), interm_size, i0.size()*i1.size()*i3.size(),
                                          1.0, data1.get(), i2.size(), transp.get(), interm_size, 0.0, interm.get(), i2.size());

              size_t iall = 0;
              for (int j013 = 0; j013 != interm_size; ++j013)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++iall)
                  interm[iall] /= min(-0.1, e0_ - (denom_->denom_xxh(j013) - eig_[j2]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i2.size(), i0.size()*i1.size()*i3.size(), interm_size,
                                          1.0, interm.get(), i2.size(), transp2.get(), interm_size, 0.0, data0.get(), i2.size());

              // sort back to the original order
              sort_indices<1,2,0,3,0,1,1,1>(data0, data1, i2.size(), i0.size(), i1.size(), i3.size());
              t->at(jst)->add_block(data1, i0, i1, i2, i3);
            }
          }
        }
      }

    } // jst loop
  } // ist loop

}

#endif
