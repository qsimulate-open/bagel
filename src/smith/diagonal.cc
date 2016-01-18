//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diagonal.cc
// Copyright (C) 2015 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <src/smith/mrci/MRCI.h>
#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/RelMRCI.h>
#include <src/smith/relcaspt2/RelCASPT2.h>
#if 0
#include <src/smith/relmrci/RelMRCI.h>
#endif

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t) const {
  for (auto& i3 : virt_) {
    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->get_size_alloc(i0, i1, i2, i3)) continue;
          unique_ptr<double[]>       data0 = t->get_block(i0, i1, i2, i3);
          const unique_ptr<double[]> data1 = t->get_block(i0, i3, i2, i1);

          sort_indices<0,3,2,1,8,1,-4,1>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                  data0[iall] *= -(eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]);
          r->add_block(data0, i0, i1, i2, i3);
        }
      }
    }
  }
}


void RelCASPT2::RelCASPT2::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t) const {
  for (auto& i3 : virt_) {
    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->get_size_alloc(i0, i1, i2, i3)) continue;
          unique_ptr<complex<double>[]> data = t->get_block(i0, i1, i2, i3);
          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                  data[iall] *= -(eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]) * 4.0;
          r->add_block(data, i0, i1, i2, i3);
        }
      }
    }
  }
}


// this function takes care of 4-external.
void MRCI::MRCI::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t) const {
  const bool diag = rdm0_->get_block()[0] == 1.0;
  if (diag) {
    for (auto& i3 : virt_) {
      for (auto& i2 : closed_) {
        for (auto& i1 : virt_) {
          for (auto& i0 : closed_) {
            // if this block is not included in the current wave function, skip it
            const size_t tsize = r->get_size_alloc(i0, i1, i2, i3);
            if (!tsize) continue;
            unique_ptr<double[]> local = r->move_block(i0, i1, i2, i3);
            unique_ptr<double[]> buf(new double[tsize]);

            for (auto& i3t : virt_) {
              for (auto& i1t : virt_) {
                unique_ptr<double[]> data0 = t->get_block(i0, i1t, i2, i3t);
                unique_ptr<double[]> data1 = t->get_block(i0, i3t, i2, i1t);
                sort_indices<0,3,2,1,8,1,-4,1>(data1, data0, i0.size(), i3t.size(), i2.size(), i1t.size());
                sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1t.size(), i2.size(), i3t.size());

                const size_t size = v2_->get_size_alloc(i1t, i1, i3t, i3);
                unique_ptr<double[]> data2 = v2_->get_block(i1t, i1, i3t, i3);
                unique_ptr<double[]> data3(new double[size]);
                sort_indices<0,2,1,3,0,1,1,1>(data2, data3, i1t.size(), i1.size(), i3t.size(), i3.size());

                dgemm_("N", "N", i0.size()*i2.size(), i1.size()*i3.size(), i1t.size()*i3t.size(),
                        1.0, data1, i0.size()*i2.size(), data3, i1t.size()*i3t.size(), 0.0, buf, i0.size()*i2.size());

                sort_indices<0,2,1,3,1,1,1,1>(buf, local, i0.size(), i2.size(), i1.size(), i3.size());
              }
            }
            r->put_block(local, i0, i1, i2, i3);
          }
        }
      }
    }
  }
}


// this function takes care of 4-external.
#if 0
void RelMRCI::RelMRCI::diagonal(shared_ptr<TATensor<std::complex<double>,4>> r, shared_ptr<const TATensor<std::complex<double>,4>> t) const {
#if 0
  const bool diag = (*rdm0_)("") == 1.0;
  if (diag)
    (*r)("c3,a4,c1,a2") += (*v2_)("a6,a2,a5,a4") * ((*t)("c1,a6,c3,a5") * 4.0);
  (*r)("c2,a3,x0,a1") += (*v2_)("a5,a1,a4,a3") * ((*t)("x1,a5,c2,a4") * ((*rdm1_)("x1,x0") * 2.0));
  (*r)("x1,a2,x0,a1") += (*v2_)("a4,a1,a3,a2") * ((*t)("x3,a3,x2,a4") * ((*rdm2_)("x3,x1,x2,x0") * 2.0));
#endif
}
#endif

#endif
