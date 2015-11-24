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
shared_ptr<TATensor<DataType,4>> SpinFreeMethod<DataType>::update_amplitude(shared_ptr<TATensor<DataType,4>> t, shared_ptr<const TATensor<DataType,4>> r) const {
  auto tt = make_shared<MultiTATensor<DataType,4>>(vector<DataType>{0.0}, vector<shared_ptr<TATensor<DataType,4>>>{t});
  // ... not good code...
  auto r0 = const_pointer_cast<TATensor<DataType,4>>(r);
  auto rr = make_shared<MultiTATensor<DataType,4>>(vector<DataType>{0.0}, vector<shared_ptr<TATensor<DataType,4>>>{r0});

  shared_ptr<MultiTATensor<DataType,4>> tnew = update_amplitude(tt, rr);
  return tnew->at(0);
}


template<typename DataType>
shared_ptr<MultiTATensor<DataType,4>> SpinFreeMethod<DataType>::update_amplitude(shared_ptr<MultiTATensor<DataType,4>> t, shared_ptr<const MultiTATensor<DataType,4>> r) const {

  auto out = t->clone();

  if (t->nref() != r->nref())
    throw logic_error("something is wrong. SpinFreeMethod::update_amplitude");

  const int fac2 = is_same<DataType,double>::value ? 1.0 : 2.0;
  const int ncore = info_->ncore() * fac2;
  const int nocc  = (info_->nclosed() + info_->nact()) * fac2;
  const int nst = t->nref();

  for (int ist = 0; ist != nst; ++ist) {
    // TODO to be fully replaced...
    auto rist = make_shared<Tensor_<DataType>>(*r->at(ist));
    auto rnew = r->at(ist);
    t->fac(ist) = 0.0;

    // not the best structure, but I am assuming that this does not take too much time...
    for (int jst = 0; jst != nst; ++jst) {
      auto tjst = make_shared<Tensor_<DataType>>(*t->at(jst));


      auto tnew = init_amplitude();
      if (ist == jst) { // AACC
        LazyTATensor<DataType,4,DenomAACC> d({virt_, virt_, closed_, closed_}, DenomAACC(eig_, nocc, nocc, ncore, ncore));
        if (is_same<DataType,double>::value)
          (*tnew)("c0,a1,c2,a3") -= ((*rnew)("c0,a1,c2,a3") * (1.0/6.0) + (*rnew)("c0,a3,c2,a1") * (1.0/12.0)) * d("a1,a3,c0,c2");
        else
          (*tnew)("c0,a1,c2,a3") -= (*rnew)("c0,a1,c2,a3") * 0.25 * d("a1,a3,c0,c2");
      }
      { // AAXX
        LazyTATensor<DataType,3,DenomGen2<1,1>> d({ortho2_, virt_, virt_}, DenomGen2<1,1>(e0_, denom_->denom_xx(), eig_, nocc, nocc));
        shared_ptr<const TATensor<DataType,3>> s = denom_->tashalf_xx({ortho2_, active_, active_});
        TATensor<DataType,3> i0(std::vector<IndexRange>{ortho2_, virt_, virt_}, true);
        TATensor<DataType,3> i1(std::vector<IndexRange>{ortho2_, virt_, virt_}, true);
        i0("o4,a1,a3") = (*rnew)("x0,a1,x2,a3") * (*s)("o4,x0,x2");
        i1("o4,a1,a3") = i0("o4,a1,a3") * d("o4,a1,a3");
        (*tnew)("x0,a1,x2,a3") -= i1("o4,a1,a3") * 0.5 * (*s)("o4,x0,x2"); // TODO conjugate s for complex cases
      }
      { // AACX
        LazyTATensor<DataType,4,DenomGen1<1,1,-1>> d({ortho1_, virt_, virt_, closed_}, DenomGen1<1,1,-1>(e0_, denom_->denom_x(), eig_, nocc, nocc, ncore));
        shared_ptr<const TATensor<DataType,2>> s = denom_->tashalf_x({ortho1_, active_});
        TATensor<DataType,4> i0(std::vector<IndexRange>{ortho1_, virt_, virt_, closed_}, true);
        TATensor<DataType,4> i1(std::vector<IndexRange>{ortho1_, virt_, virt_, closed_}, true);
        if (is_same<DataType,double>::value)
          i0("o4,a1,a3,c2") = ((*rnew)("c2,a3,x0,a1")*(2.0/3.0) + (*rnew)("c2,a1,x0,a3")*(1.0/3.0)) * (*s)("o4,x0");
        else
          i0("o4,a1,a3,c2") = (*rnew)("c2,a3,x0,a1") * 0.5 * (*s)("o4,x0");
        i1("o4,a1,a3,c2") = i0("o4,a1,a3,c2") * d("o4,a1,a3,c2");
        (*tnew)("x0,a1,c2,a3") -= i1("o4,a1,a3,c2") * (*s)("o4,x0"); // TODO conjugate s for complex cases
      }
      { // AXCC
        LazyTATensor<DataType,4,DenomGen1<1,-1,-1>> d({ortho1_, virt_, closed_, closed_}, DenomGen1<1,-1,-1>(e0_, denom_->denom_h(), eig_, nocc, ncore, ncore));
        shared_ptr<const TATensor<DataType,2>> s = denom_->tashalf_h({ortho1_, active_});
        TATensor<DataType,4> i0(std::vector<IndexRange>{ortho1_, virt_, closed_, closed_}, true);
        TATensor<DataType,4> i1(std::vector<IndexRange>{ortho1_, virt_, closed_, closed_}, true);
        if (is_same<DataType,double>::value)
          i0("o4,a1,c0,c2") = ((*rnew)("c2,x3,c0,a1")*(2.0/3.0) + (*rnew)("c0,x3,c2,a1")*(1.0/3.0)) * (*s)("o4,x3");
        else
          i0("o4,a1,c0,c2") = (*rnew)("c2,x3,c0,a1") * 0.5 * (*s)("o4,x3");
        i1("o4,a1,c0,c2") = i0("o4,a1,c0,c2") * d("o4,a1,c0,c2");
        (*tnew)("c0,a1,c2,x3") -= i1("o4,a1,c0,c2") * (*s)("o4,x3"); // TODO conjugate s for complex cases
      }
      { // XXCC
        LazyTATensor<DataType,3,DenomGen2<-1,-1>> d({ortho2_, closed_, closed_}, DenomGen2<-1,-1>(e0_, denom_->denom_hh(), eig_, ncore, ncore));
        shared_ptr<const TATensor<DataType,3>> s = denom_->tashalf_hh({ortho2_, active_, active_});
        TATensor<DataType,3> i0(std::vector<IndexRange>{ortho2_, closed_, closed_}, true);
        TATensor<DataType,3> i1(std::vector<IndexRange>{ortho2_, closed_, closed_}, true);
        i0("o4,c0,c2") = (*rnew)("c0,x1,c2,x3") * (*s)("o4,x1,x3");
        i1("o4,c0,c2") = i0("o4,c0,c2") * d("o4,c0,c2");
        (*tnew)("c0,x1,c2,x3") -= i1("o4,c0,c2") * 0.5 * (*s)("o4,x1,x3"); // TODO conjugate s for complex cases
      }
      { // XXCX
        LazyTATensor<DataType,2,DenomGen3<-1>> d({ortho3_, closed_}, DenomGen3<-1>(e0_, denom_->denom_xxh(), eig_, ncore));
        shared_ptr<const TATensor<DataType,4>> s = denom_->tashalf_xxh({ortho3_, active_, active_, active_});
        TATensor<DataType,2> i0(std::vector<IndexRange>{ortho3_, closed_}, true);
        TATensor<DataType,2> i1(std::vector<IndexRange>{ortho3_, closed_}, true);
        i0("o4,c2") = (*rnew)("c2,x3,x0,x1") * (*s)("o4,x0,x1,x3");
        i1("o4,c2") = i0("o4,c2") * d("o4,c2");
        (*tnew)("x0,x1,c2,x3") -= i1("o4,c2") * (*s)("o4,x0,x1,x3"); // TODO conjugate s for complex cases
      }
      { // AXXX
        LazyTATensor<DataType,2,DenomGen3<1>> d({ortho3_, virt_}, DenomGen3<1>(e0_, denom_->denom_xhh(), eig_, nocc));
        shared_ptr<const TATensor<DataType,4>> s = denom_->tashalf_xhh({ortho3_, active_, active_, active_});
        TATensor<DataType,2> i0(std::vector<IndexRange>{ortho3_, virt_}, true);
        TATensor<DataType,2> i1(std::vector<IndexRange>{ortho3_, virt_}, true);
        i0("o4,a1") = (*rnew)("x2,x3,x0,a1") * (*s)("o4,x0,x2,x3");
        i1("o4,a1") = i0("o4,a1") * d("o4,a1");
        (*tnew)("x0,a1,x2,x3") -= i1("o4,a1") * (*s)("o4,x0,x2,x3"); // TODO conjugate s for complex cases
      }

      for (auto& i3 : active_) {
      for (auto& i2 : active_) {
      if (is_same<DataType,double>::value) {
        assert(denom_->shalf_xh());
        const size_t interm_size = denom_->shalf_xh()->ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&nclo,&nact,&interm_size, this](const int i, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*interm_size*2]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k) {
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclo)+(j3-nclo)*nact + 2*i*nact*nact),
                     interm_size, out.get()+interm_size*k);
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclo)+(j3-nclo)*nact + (2*i+1)*nact*nact),
                     interm_size, out.get()+interm_size*(k+I2.size()*I3.size()));
            }
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ist, i2, i3);

        for (auto& i3t : active_) {
        for (auto& i2t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jst, i2t, i3t);
          blas::conj_n(transp2.get(), i2t.size()*i3t.size()*interm_size*2);

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = rist->get_size_alloc(i2, i3, i0, i1);
              const size_t blocksizet = tjst->get_size_alloc(i0, i1, i2t, i3t);
//            const size_t blocksizet = rjst->get_size_alloc(i2t, i3t, i0, i1);
              if (!blocksize || !blocksizet) continue;
              assert(blocksize == rist->get_size_alloc(i0, i3, i2, i1));
              unique_ptr<DataType[]> data0 = rist->get_block(i2, i3, i0, i1);
              unique_ptr<DataType[]> data1 = rist->get_block(i0, i3, i2, i1);

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
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom_xh(j23) + eig_[j1] - eig_[j0]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2t.size()*i3t.size()*2, interm_size,
                                          1.0, interm.get(), i0.size()*i1.size(), transp2.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size());

              // sort back to the original order
              unique_ptr<DataType[]> data3(new DataType[blocksizet]);
              unique_ptr<DataType[]> data4(new DataType[blocksizet]);
              copy_n(data2.get(), blocksizet, data3.get());
              sort_indices<2,1,0,3,0,1,1,1>(data2.get()+blocksizet, data4.get(), i0.size(), i1.size(), i2t.size(), i3t.size());
              tjst->add_block(data3, i0, i1, i2t, i3t);
              tjst->add_block(data4, i2t, i1, i0, i3t);
            }
          }
        }
        }
      } else {
        assert(denom_->shalf_xh());
        const size_t interm_size = denom_->shalf_xh()->ndim();
        const int nact = info_->nact() * fac2;
        const int nclo = info_->nclosed() * fac2;
        auto create_transp = [&,this](const int i, const Index& I2, const Index& I3) {
          unique_ptr<DataType[]> out(new DataType[I2.size()*I3.size()*interm_size]);
          for (int j3 = I3.offset(), k = 0; j3 != I3.offset()+I3.size(); ++j3)
            for (int j2 = I2.offset(); j2 != I2.offset()+I2.size(); ++j2, ++k)
              copy_n(denom_->shalf_xh()->element_ptr(0, (j2-nclo)+(j3-nclo)*nact + i*nact*nact),
                     interm_size, out.get()+interm_size*k);
          return move(out);
        };
        unique_ptr<DataType[]> transp = create_transp(ist, i2, i3);

        for (auto& i3t : active_) {
        for (auto& i2t : active_) {
          unique_ptr<DataType[]> transp2 = create_transp(jst, i2t, i3t);
          blas::conj_n(transp2.get(), i2t.size()*i3t.size()*interm_size);

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = rist->get_size_alloc(i2, i3, i0, i1);
              const size_t blocksizet = tjst->get_size_alloc(i0, i1, i2t, i3t);
//            const size_t blocksizet = rjst->get_size_alloc(i2t, i3t, i0, i1);
              if (!blocksize || !blocksizet) continue;
              unique_ptr<DataType[]> data0 = rist->get_block(i2, i3, i0, i1);

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
                    interm[iall] /= min(-0.1, e0_ - (denom_->denom_xh(j23) + eig_[j1] - eig_[j0]));

              // move back to non-orthogonal basis
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, i0.size()*i1.size(), i2t.size()*i3t.size(), interm_size,
                                          1.0, interm.get(), i0.size()*i1.size(), transp2.get(), interm_size, 0.0, data2.get(), i0.size()*i1.size());

              // sort back to the original order
              unique_ptr<DataType[]> data3(new DataType[blocksizet]);
              copy_n(data2.get(), blocksizet, data3.get());
              tjst->add_block(data3, i0, i1, i2t, i3t);
            }
          }
        }
        }
      }
      }
      }

      auto tmp = tjst->template tiledarray<4>();
      (*tmp)("c0,a1,c2,a3") += (*tnew)("c0,a1,c2,a3");
      (*tmp)("x0,a1,x2,a3") += (*tnew)("x0,a1,x2,a3");
      (*tmp)("x0,a1,c2,a3") += (*tnew)("x0,a1,c2,a3");
      (*tmp)("c0,a1,c2,x3") += (*tnew)("c0,a1,c2,x3");
      (*tmp)("c0,x1,c2,x3") += (*tnew)("c0,x1,c2,x3");
      (*tmp)("x0,x1,c2,x3") += (*tnew)("x0,x1,c2,x3");
      (*tmp)("x0,a1,x2,x3") += (*tnew)("x0,a1,x2,x3");
      (*out)[jst] = tmp;
    } // jst loop
  } // ist loop

  return out;
}

#endif
