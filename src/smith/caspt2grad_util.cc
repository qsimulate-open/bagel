//
// BAGEL - Parallel electron correlation program.
// Filename: caspt2grad_util.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <bagel_config.h>
#include <src/smith/caspt2grad.h>
#ifdef COMPILE_SMITH
#include <ga.h>
#endif

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<DFFullDist> CASPT2Grad::contract_D1(shared_ptr<const DFFullDist> full) const {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nocc = ref_->nocc();
  const int nall = nocc + ref_->nvirt();
  const int n = nocc-ncore_;

  shared_ptr<DFFullDist> out = full->clone();
  double* optr = out->block(0)->data();
  const double* iptr = full->block(0)->data();
  const size_t asize = full->block(0)->asize();

  // first create a tensor for full.
  // Aux index blocking
  shared_ptr<const StaticDist> adist = ref_->geom()->df()->adist_now();
  const IndexRange aux(adist);

  vector<IndexRange> ind = d2_->indexrange();
  assert(ind[0] == ind[2] && ind[1] == ind[3]);

  // create a GA array
  auto input = make_shared<Tensor>(vector<IndexRange>{aux, ind[0], ind[1]});
  input->allocate();
  auto output = input->clone();

  // fill into input
  for (auto& a : aux) {
    tuple<size_t, size_t> loc = adist->locate(a.offset());
    if (get<0>(loc) != mpi__->rank())
      continue;
    for (auto& j : ind[1])
      for (auto& i : ind[0]) {
        unique_ptr<double[]> buf(new double[a.size()*i.size()*j.size()]);
        for (int ej = 0; ej != j.size(); ++ej)
          for (int ei = 0; ei != i.size(); ++ei)
            copy_n(iptr + get<1>(loc)+asize*(i.offset()+ei+full->nocc1()*(j.offset()+ej)), a.size(), buf.get()+a.size()*(ei+i.size()*ej));
        input->put_block(buf, a, i, j);
      }
  }
  GA_Sync();

  // next contract
  for (auto& j : ind[1])
    for (auto& i : ind[0])
      for (auto& a : aux) {
        if (!output->is_local(a, i, j)) continue;
        unique_ptr<double[]> buf(new double[a.size()*i.size()*j.size()]);
        fill_n(buf.get(), a.size()*i.size()*j.size(), 0.0);
        for (auto& l : ind[1]) {
          for (auto& k : ind[0]) {
            unique_ptr<double[]> in = input->get_block(a, k, l);
            if (d2_->exists(k, l, i, j)) {
              unique_ptr<double[]> d = d2_->get_block(k, l, i, j);
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasNoTrans, a.size(), i.size()*j.size(), k.size()*l.size(),
                                          1.0, in.get(), a.size(), d.get(), k.size()*l.size(), 1.0, buf.get(), a.size());
            } else if (d2_->exists(i, j, k, l)) {
              unique_ptr<double[]> d = d2_->get_block(i, j, k, l);
              btas::gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, a.size(), i.size()*j.size(), k.size()*l.size(),
                                          1.0, in.get(), a.size(), d.get(), i.size()*j.size(), 1.0, buf.get(), a.size());
            }
          }
        }
        output->add_block(buf, a, i, j);
      }
  GA_Sync();

  // finally fill back into output
  for (auto& a : aux) {
    tuple<size_t, size_t> loc = adist->locate(a.offset());
    if (get<0>(loc) != mpi__->rank())
      continue;
    for (auto& j : ind[1])
      for (auto& i : ind[0]) {
        unique_ptr<double[]> buf = output->get_block(a, i, j);
        for (int ej = 0; ej != j.size(); ++ej)
          for (int ei = 0; ei != i.size(); ++ei)
            copy_n(buf.get()+a.size()*(ei+i.size()*ej), a.size(), optr + get<1>(loc)+asize*(i.offset()+ei+full->nocc1()*(j.offset()+ej)));
      }
  }
  GA_Sync();

  // one body part
  {
    auto is_cc   = [&](const int& i) { return i < nclosed; };
    auto is_act  = [&](const int& i) { return i >= nclosed && i < nocc; };
    auto is_virt = [&](const int& i) { return i >= nocc; };
    auto is_exc  = [&](const int& i, const int& j) { return is_virt(i) || (is_act(i) && is_cc(j)); };

    for (int t = 0; t != nall; ++t)
      for (int l = 0; l != nocc; ++l)
        // c(cc)x, c(cc)a, x(cc)a
        if (is_exc(t, l))
          for (int s = 0; s != nclosed; ++s) {
            blas::ax_plus_y_n(2.0*d11_->element(l,t), iptr + asize*(l+nocc*t), asize, optr + asize*(s+nocc*s));
            blas::ax_plus_y_n(   -d11_->element(l,t), iptr + asize*(s+nocc*t), asize, optr + asize*(l+nocc*s));
            blas::ax_plus_y_n(2.0*d11_->element(l,t), iptr + asize*(s+nocc*s), asize, optr + asize*(l+nocc*t));
            blas::ax_plus_y_n(   -d11_->element(l,t), iptr + asize*(l+nocc*s), asize, optr + asize*(s+nocc*t));
          }
  }
  return out;
#else
  return full->clone(); // dummy
#endif
}
