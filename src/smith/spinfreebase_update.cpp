//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spinfreebase_denom.cpp
// Copyright (C) 2015 Toru Shiozaki
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

    if (t->at(ist)) {
      const double e0loc = e0all_[ist] - e0_;
      for (auto& i3 : virt_) {
        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              if (!t->at(ist)->is_local(i0, i1, i2, i3) || !r->at(ist)->get_size(i0, i1, i2, i3)) continue;
              unique_ptr<DataType[]>       data0 = r->at(ist)->get_block(i0, i1, i2, i3);

              // this is an inverse of the overlap.
              if (is_same<DataType,double>::value) {
                const unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i0, i3, i2, i1);
                sort_indices<0,3,2,1,2,12,1,12>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
              } else {
                blas::scale_n(0.25, data0.get(), r->at(ist)->get_size(i0, i1, i2, i3));
              }
              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                      // note that e0 is cancelled by another term
                      data0[iall] /= eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1] - e0loc;
              t->at(ist)->add_block(data0, i0, i1, i2, i3);
            }
          }
        }
      }
    }

    // not the best structure, but I am assuming that this does not take too much time...
    for (int jst = 0; jst != nst; ++jst) {
      if (!t->at(jst) || !r->at(ist)) continue;

      {
      const ViewType ishalf = denom_->shalf("xx", ist);
      const ViewType jshalf = denom_->shalf("xx", jst);
      for (auto& i2 : active_) {
      for (auto& i0 : active_) {
        // trans is the transformation matrix
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I0, const Index& I2) {
          unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*interm_size]);
          for (int j2 = I2.offset(), k = 0; j2 != I2.offset()+I2.size(); ++j2)
            for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
              copy_n(shalf.element_ptr(0, (j0-nclo)+(j2-nclo)*nact), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i0, i2);

        for (auto& i2t : active_) {
        for (auto& i0t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i0t, i2t);

          for (auto& i3 : virt_) {
            for (auto& i1 : virt_) {
              if (!t->at(jst)->is_local(i0t, i1, i2t, i3)) continue;
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size(i0, i1, i2, i3);
              const size_t blocksizet = r->at(jst)->get_size(i0t, i1, i2t, i3);
              if (!blocksize || !blocksizet) continue;
              // data0 is the source area
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i0, i1, i2, i3);
              unique_ptr<DataType[]> data1(new DataType[max(blocksize,blocksizet)]);
              // sort. Active indices run faster
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i1.size()*i3.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, interm_size, i1.size()*i3.size(), i0.size()*i2.size(),
                                          1.0, transp.get(), interm_size, data1.get(), i0.size()*i2.size(), 0.0, interm.get(), interm_size);

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j02 = 0; j02 != interm_size; ++j02, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom("xx", ist, j02) + eig_[j3] + eig_[j1]));

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              // TODO check for complex cases
              btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0t.size()*i2t.size(), i1.size()*i3.size(), interm_size,
                                          0.5, transp2.get(), interm_size, interm.get(), interm_size, 0.0, data1.get(), i0t.size()*i2t.size());

              // sort back to the original order
              unique_ptr<DataType[]> data2(new DataType[blocksizet]);
              sort_indices<0,2,1,3,0,1,1,1>(data1, data2, i0t.size(), i2t.size(), i1.size(), i3.size());
              t->at(jst)->add_block(data2, i0t, i1, i2t, i3);
            }
          }
        }
        }
      }
      }
      }

      {
      const ViewType ishalf = denom_->shalf("x", ist);
      const ViewType jshalf = denom_->shalf("x", jst);
      for (auto& i0 : active_) {
        // trans is the transformation matrix
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I0) {
          unique_ptr<DataType[]> out(new DataType[I0.size()*interm_size]);
          for (int j0 = I0.offset(), k = 0; j0 != I0.offset()+I0.size(); ++j0, ++k)
            copy_n(shalf.element_ptr(0,j0-nclo), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i0);

        for (auto& i0t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i0t);

          for (auto& i3 : virt_) {
            for (auto& i2 : closed_) {
              for (auto& i1 : virt_) {
                if (!t->at(jst)->is_local(i0t, i1, i2, i3)) continue;
                const size_t blocksize = r->at(ist)->get_size(i2, i3, i0, i1);
                const size_t blocksizet = r->at(jst)->get_size(i2, i3, i0t, i1);
                if (!blocksize || !blocksizet) continue;

                unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
                unique_ptr<DataType[]> data2(new DataType[blocksize]);
                sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
                if (is_same<DataType,double>::value) {
                  assert(r->at(ist)->get_size(i2, i1, i0, i3));
                  const unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i2, i1, i0, i3);
                  sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());
                } else {
                  blas::scale_n(0.5, data2.get(), blocksize);
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
                        interm[iall] /= min(-0.1, e0_ - (denom_->denom("x", ist, j0) + eig_[j3] - eig_[j2] + eig_[j1]));

                // move back to non-orthogonal basis
                unique_ptr<DataType[]> data3(new DataType[blocksizet]);
                btas::gemm_impl<true>::call(CblasColMajor, CblasConjTrans, CblasNoTrans, i0t.size(), i1.size()*i2.size()*i3.size(), interm_size,
                                            1.0, transp2.get(), interm_size, interm.get(), interm_size, 0.0, data3.get(), i0t.size());

                t->at(jst)->add_block(data3, i0t, i1, i2, i3);
              }
            }
          }
        }
      }
      }

      {
      const ViewType ishalf = denom_->shalf("h", ist);
      const ViewType jshalf = denom_->shalf("h", jst);
      for (auto& i3 : active_) {
        // trans is the transformation matrix
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3, ++k)
            copy_n(shalf.element_ptr(0,j3-nclo), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i3);

        for (auto& i3t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i3t);
          blas::conj_n(transp2.get(), i3t.size()*interm_size);

          for (auto& i2 : closed_) {
            for (auto& i1 : virt_) {
              for (auto& i0 : closed_) {
                if (!t->at(jst)->is_local(i0, i1, i2, i3t)) continue;
                const size_t blocksize = r->at(ist)->get_size(i2, i3, i0, i1);
                const size_t blocksizet = r->at(jst)->get_size(i2, i3t, i0, i1);
                if (!blocksize || !blocksizet) continue;

                assert(r->at(ist)->get_size(i0, i3, i2, i1));
                unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
                unique_ptr<DataType[]> data2(new DataType[blocksize]);
                sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
                if (is_same<DataType,double>::value) {
                  const unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i0, i3, i2, i1);
                  sort_indices<0,3,2,1,2,3,1,3>(data1, data2, i0.size(), i3.size(), i2.size(), i1.size());
                } else {
                  blas::scale_n(0.5, data2.get(), blocksize);
                }

                // move to orthogonal basis
                unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*i2.size()*interm_size]);
                btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size()*i2.size(), interm_size, i3.size(),
                                            1.0, data2.get(), i0.size()*i1.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size()*i2.size());

                size_t iall = 0;
                for (int j3 = 0; j3 != interm_size; ++j3)
                  for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                        interm[iall] /= min(-0.1, e0_ - (denom_->denom("h", ist, j3) - eig_[j2] + eig_[j1] - eig_[j0]));

                // move back to non-orthogonal basis
                unique_ptr<DataType[]> data3(new DataType[blocksizet]);
                btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size()*i2.size(), i3t.size(), interm_size,
                                            1.0, interm.get(), i0.size()*i1.size()*i2.size(), transp2.get(), interm_size, 0.0, data3.get(), i0.size()*i1.size()*i2.size());

                t->at(jst)->add_block(data3, i0, i1, i2, i3t);
              }
            }
          }
        }
      }
      }

      {
      const ViewType ishalf = denom_->shalf("hh", ist);
      const ViewType jshalf = denom_->shalf("hh", jst);
      for (auto& i3 : active_) {
      for (auto& i1 : active_) {
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I1, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I1.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1, ++k)
              copy_n(shalf.element_ptr(0,(j1-nclo)+(j3-nclo)*nact), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i1, i3);

        for (auto& i3t : active_) {
        for (auto& i1t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i1t, i3t);
          blas::conj_n(transp2.get(), i1t.size()*i3t.size()*interm_size);

          for (auto& i2 : closed_) {
            for (auto& i0 : closed_) {
              if (!t->at(jst)->is_local(i0, i1t, i2, i3t)) continue;
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size(i0, i1, i2, i3);
              const size_t blocksizet = r->at(jst)->get_size(i0, i1t, i2, i3t);
              if (!blocksize || !blocksizet) continue;
              // data0 is the source area
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i0, i1, i2, i3);
              unique_ptr<DataType[]> data1(new DataType[max(blocksize,blocksizet)]);
              // sort. Active indices run slower
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i0.size()*i2.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i2.size(), interm_size, i1.size()*i3.size(),
                                          1.0, data1.get(), i0.size()*i2.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i2.size());

              size_t iall = 0;
              for (int j13 = 0; j13 != interm_size; ++j13)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom("hh", ist, j13) - eig_[j2] - eig_[j0]));

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i2.size(), i1t.size()*i3t.size(), interm_size,
                                          0.5, interm.get(), i0.size()*i2.size(), transp2.get(), interm_size, 0.0, data1.get(), i0.size()*i2.size());

              // sort back to the original order
              unique_ptr<DataType[]> data2(new DataType[blocksizet]);
              sort_indices<0,2,1,3,0,1,1,1>(data1, data2, i0.size(), i2.size(), i1t.size(), i3t.size());
              t->at(jst)->add_block(data2, i0, i1t, i2, i3t);
            }
          }
        }
        }
      }
      }
      }

      {
      const ViewType ishalf = denom_->shalf("xh", ist);
      const ViewType jshalf = denom_->shalf("xh", jst);
      for (auto& i3 : active_) {
      for (auto& i2 : active_) {
      if (is_same<DataType,double>::value) {
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*interm_size*2]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
              copy_n(shalf.element_ptr(0, (j2-nclo)+(j3-nclo)*nact), interm_size, out.get()+interm_size*k);
              copy_n(shalf.element_ptr(0, (j2-nclo)+(j3-nclo)*nact + nact*nact), interm_size, out.get()+interm_size*(k+I2.size()*I3.size()));
            }
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i2, i3);

        for (auto& i3t : active_) {
        for (auto& i2t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i2t, i3t);
          blas::conj_n(transp2.get(), i2t.size()*i3t.size()*interm_size*2);

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              if (!t->at(jst)->is_local(i0, i1, i2t, i3t)) continue;
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size(i2, i3, i0, i1);
              const size_t blocksizet = r->at(jst)->get_size(i2t, i3t, i0, i1);
              if (!blocksize || !blocksizet) continue;
              assert(blocksize == r->at(ist)->get_size(i0, i3, i2, i1));
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
              unique_ptr<DataType[]> data1 = r->at(ist)->get_block(i0, i3, i2, i1);

              unique_ptr<DataType[]> data2(new DataType[max(blocksize,blocksizet)*2]);
              // sort. Active indices run slower
              sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
              sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size(), interm_size, i2.size()*i3.size()*2,
                                          1.0, data2.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size());

              size_t iall = 0;
              for (int j23 = 0; j23 != interm_size; ++j23)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom("xh", ist, j23) + eig_[j1] - eig_[j0]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2t.size()*i3t.size()*2, interm_size,
                                          1.0, interm.get(), i0.size()*i1.size(), transp2.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size());

              // sort back to the original order
              unique_ptr<DataType[]> data3(new DataType[blocksizet]);
              unique_ptr<DataType[]> data4(new DataType[blocksizet]);
              copy_n(data2.get(), blocksizet, data3.get());
              sort_indices<2,1,0,3,0,1,1,1>(data2.get()+blocksizet, data4.get(), i0.size(), i1.size(), i2t.size(), i3t.size());
              t->at(jst)->add_block(data3, i0, i1, i2t, i3t);
              t->at(jst)->add_block(data4, i2t, i1, i0, i3t);
            }
          }
        }
        }
      } else {
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&,this](const ViewType shalf, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k)
              copy_n(shalf.element_ptr(0, (j2-nclo)+(j3-nclo)*nact), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i2, i3);

        for (auto& i3t : active_) {
        for (auto& i2t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i2t, i3t);
          blas::conj_n(transp2.get(), i2t.size()*i3t.size()*interm_size);

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              if (!t->at(jst)->is_local(i0, i1, i2t, i3t)) continue;
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->at(ist)->get_size(i2, i3, i0, i1);
              const size_t blocksizet = r->at(jst)->get_size(i2t, i3t, i0, i1);
              if (!blocksize || !blocksizet) continue;
              unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);

              unique_ptr<DataType[]> data2(new DataType[max(blocksize,blocksizet)]);
              // sort. Active indices run slower
              sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get(), i2.size(), i3.size(), i0.size(), i1.size());
              // intermediate area
              unique_ptr<DataType[]> interm(new DataType[i0.size()*i1.size()*interm_size]);

              // move to orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i0.size()*i1.size(), interm_size, i2.size()*i3.size(),
                                          1.0, data2.get(), i0.size()*i1.size(), transp.get(), interm_size, 0.0, interm.get(), i0.size()*i1.size());

              size_t iall = 0;
              for (int j23 = 0; j23 != interm_size; ++j23)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom("xh", ist, j23) + eig_[j1] - eig_[j0]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2t.size()*i3t.size(), interm_size,
                                          1.0, interm.get(), i0.size()*i1.size(), transp2.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size());

              // sort back to the original order
              unique_ptr<DataType[]> data3(new DataType[blocksizet]);
              copy_n(data2.get(), blocksizet, data3.get());
              t->at(jst)->add_block(data3, i0, i1, i2t, i3t);
            }
          }
        }
        }
      }
      }
      }
      }

      {
      const ViewType ishalf = denom_->shalf("xxh", ist);
      const ViewType jshalf = denom_->shalf("xxh", jst);
      for (auto& i3 : active_) {
      for (auto& i2 : active_) {
      for (auto& i0 : active_) {
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I0, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I0.size()*I2.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2)
              for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                copy_n(shalf.element_ptr(0,j0-nclo+nact*(j2-nclo+nact*(j3-nclo))), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i0, i2, i3);

        for (auto& i3t : active_) {
        for (auto& i2t : active_) {
        for (auto& i0t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i0t, i2t, i3t);
          blas::conj_n(transp2.get(), i0t.size()*i2t.size()*i3t.size()*interm_size);

          for (auto& i1 : virt_) {
            if (!t->at(jst)->is_local(i0t, i1, i2t, i3t)) continue;
            // if this block is not included in the current wave function, skip it
            const size_t blocksize = r->at(ist)->get_size(i2, i3, i0, i1);
            const size_t blocksizet = r->at(jst)->get_size(i2t, i3t, i0t, i1);
            if (!blocksize || !blocksizet) continue;
            // data0 is the source area
            unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data1(new DataType[max(blocksize,blocksizet)]);
            // sort. Active indices run slower
            sort_indices<3,2,0,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
            // intermediate area
            unique_ptr<DataType[]> interm(new DataType[i1.size()*interm_size]);

            // move to orthogonal basis
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i1.size(), interm_size, i0.size()*i2.size()*i3.size(),
                                        1.0, data1.get(), i1.size(), transp.get(), interm_size, 0.0, interm.get(), i1.size());

            size_t iall = 0;
            for (int j123 = 0; j123 != interm_size; ++j123)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++iall)
                interm[iall] /= min(-0.1, e0_ - (denom_->denom("xxh", ist, j123) + eig_[j1]));

            // move back to non-orthogonal basis
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i1.size(), i0t.size()*i2t.size()*i3t.size(), interm_size,
                                        1.0, interm.get(), i1.size(), transp2.get(), interm_size, 0.0, data1.get(), i1.size());

            // sort back to the original order
            unique_ptr<DataType[]> data2(new DataType[blocksizet]);
            sort_indices<1,0,2,3,0,1,1,1>(data1, data2, i1.size(), i0t.size(), i2t.size(), i3t.size());
            t->at(jst)->add_block(data2, i0t, i1, i2t, i3t);
          }
        }
        }
        }
      }
      }
      }
      }

      {
      const ViewType ishalf = denom_->shalf("xhh", ist);
      const ViewType jshalf = denom_->shalf("xhh", jst);
      for (auto& i3 : active_) {
      for (auto& i1 : active_) {
      for (auto& i0 : active_) {
        const size_t interm_size = ishalf.ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const ViewType shalf, const Index& I0, const Index& I1, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I0.size()*I1.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j1 = I1.offset(); j1 != I1.offset()+I1.size(); ++j1)
              for (int j0 = I0.offset(); j0 != I0.offset()+I0.size(); ++j0, ++k)
                copy_n(shalf.element_ptr(0,j0-nclo+nact*(j1-nclo+nact*(j3-nclo))), interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ishalf, i0, i1, i3);

        for (auto& i3t : active_) {
        for (auto& i1t : active_) {
        for (auto& i0t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jshalf, i0t, i1t, i3t);
          blas::conj_n(transp2.get(), i0t.size()*i1t.size()*i3t.size()*interm_size);

          for (auto& i2 : closed_) {
            if (!t->at(jst)->is_local(i0t, i1t, i2, i3t)) continue;
            // if this block is not included in the current wave function, skip it
            const size_t blocksize = r->at(ist)->get_size(i2, i3, i0, i1);
            const size_t blocksizet = r->at(jst)->get_size(i2, i3t, i0t, i1t);
            if (!blocksize || !blocksizet) continue;
            // data0 is the source area
            unique_ptr<DataType[]> data0 = r->at(ist)->get_block(i2, i3, i0, i1);
            unique_ptr<DataType[]> data1(new DataType[max(blocksize,blocksizet)]);
            // sort. Active indices run slower
            sort_indices<0,2,3,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
            // intermediate area
            unique_ptr<DataType[]> interm(new DataType[i2.size()*interm_size]);

            // move to orthogonal basis
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, i2.size(), interm_size, i0.size()*i1.size()*i3.size(),
                                        1.0, data1.get(), i2.size(), transp.get(), interm_size, 0.0, interm.get(), i2.size());

            size_t iall = 0;
            for (int j013 = 0; j013 != interm_size; ++j013)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++iall)
                interm[iall] /= min(-0.1, e0_ - (denom_->denom("xhh", ist, j013) - eig_[j2]));

            // move back to non-orthogonal basis
            btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i2.size(), i0t.size()*i1t.size()*i3t.size(), interm_size,
                                        1.0, interm.get(), i2.size(), transp2.get(), interm_size, 0.0, data1.get(), i2.size());

            // sort back to the original order
            unique_ptr<DataType[]> data2(new DataType[blocksizet]);
            sort_indices<1,2,0,3,0,1,1,1>(data1, data2, i2.size(), i0t.size(), i1t.size(), i3t.size());
            t->at(jst)->add_block(data2, i0t, i1t, i2, i3t);
          }
        }
        }
        }
      }
      }
      }
      }

    } // jst loop
  } // ist loop

  mpi__->barrier();
}

#endif
