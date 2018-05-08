//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diagonal.cc
// Copyright (C) 2015 Toru Shiozaki
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
#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relcaspt2/RelCASPT2.h>
#include <src/smith/casa/CASA.h>
#include <src/smith/relcasa/RelCASA.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t, const bool diag) const {
  double sum = 0.0;
  for (auto& i1 : active_) {
    for (auto& i0 : active_) {
      if (f1_->is_local(i0, i1)) {
        unique_ptr<double[]> fdata = f1_->get_block(i0, i1);
        unique_ptr<double[]> rdata = rdm1_->get_block(i0, i1);
        sum += blas::dot_product_noconj(fdata.get(), i0.size()*i1.size(), rdata.get());
      }
    }
  }
  mpi__->allreduce(&sum, 1);

  const double e0loc = sum - (diag ? e0_ : 0.0);
  for (auto& i3 : virt_) {
    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->is_local(i0, i1, i2, i3) || !r->get_size(i0, i1, i2, i3)) continue;
          unique_ptr<double[]>       data0 = t->get_block(i0, i1, i2, i3);
          const unique_ptr<double[]> data1 = t->get_block(i0, i3, i2, i1);

          sort_indices<0,3,2,1,8,1,-4,1>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
          if (diag) {
            size_t iall = 0;
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    data0[iall] *= e0loc - (eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]);
          } else {
            blas::scale_n(e0loc, data0.get(), i0.size()*i3.size()*i2.size()*i1.size());
          }
          r->add_block(data0, i0, i1, i2, i3);
        }
      }
    }
  }
  mpi__->barrier();
}


void RelCASPT2::RelCASPT2::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t, const bool diag) const {
  complex<double> sum = 0.0;
  for (auto& i1 : active_) {
    for (auto& i0 : active_) {
      if (f1_->is_local(i0, i1)) {
        unique_ptr<complex<double>[]> fdata = f1_->get_block(i0, i1);
        unique_ptr<complex<double>[]> rdata = rdm1_->get_block(i0, i1);
        sum += blas::dot_product_noconj(fdata.get(), i0.size()*i1.size(), rdata.get());
      }
    }
  }
  mpi__->allreduce(&sum, 1);

  const complex<double> e0loc = sum - (diag ? e0_ : 0.0);
  for (auto& i3 : virt_) {
    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->is_local(i0, i1, i2, i3) || !r->get_size(i0, i1, i2, i3)) continue;
          unique_ptr<complex<double>[]> data = t->get_block(i0, i1, i2, i3);
          if (diag) {
            size_t iall = 0;
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    data[iall] *= (e0loc -(eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1])) * 4.0;
          } else {
            blas::scale_n(4.0*e0loc, data.get(), i0.size()*i3.size()*i2.size()*i1.size());
          }
          r->add_block(data, i0, i1, i2, i3);
        }
      }
    }
  }
  mpi__->barrier();
}


void CASA::CASA::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t, const bool diag) const {
  // Not currently using optimized code for diagonal contributions
}


void RelCASA::RelCASA::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t, const bool diag) const {
  // Not currently using optimized code for diagonal contributions
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
            if (!r->is_local(i0, i1, i2, i3) || !r->get_size(i0, i1, i2, i3)) continue;
            const size_t tsize = r->get_size(i0, i1, i2, i3);
            unique_ptr<double[]> local(new double[tsize]);
            unique_ptr<double[]> buf(new double[tsize]);
            fill_n(buf.get(), tsize, 0.0);

            for (auto& i3t : virt_) {
              for (auto& i1t : virt_) {
                unique_ptr<double[]> data0 = t->get_block(i0, i1t, i2, i3t);
                unique_ptr<double[]> data1 = t->get_block(i0, i3t, i2, i1t);
                sort_indices<0,3,2,1,8,1,-4,1>(data1, data0, i0.size(), i3t.size(), i2.size(), i1t.size());
                sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1t.size(), i2.size(), i3t.size());

                unique_ptr<double[]> data2 = v2_->get_block(i1t, i1, i3t, i3);
                unique_ptr<double[]> data3(new double[v2_->get_size(i1t, i1, i3t, i3)]);
                sort_indices<0,2,1,3,0,1,1,1>(data2, data3, i1t.size(), i1.size(), i3t.size(), i3.size());

                dgemm_("N", "N", i0.size()*i2.size(), i1.size()*i3.size(), i1t.size()*i3t.size(),
                        1.0, data1, i0.size()*i2.size(), data3, i1t.size()*i3t.size(), 1.0, buf, i0.size()*i2.size());
              }
            }
            sort_indices<0,2,1,3,0,1,1,1>(buf, local, i0.size(), i2.size(), i1.size(), i3.size());
            r->add_block(local, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  mpi__->barrier();
}


// this function takes care of 4-external.
void RelMRCI::RelMRCI::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t) const {
  const bool diag = rdm0_->get_block()[0] == 1.0;
  if (diag) {
    for (auto& i3 : virt_) {
      for (auto& i2 : closed_) {
        for (auto& i1 : virt_) {
          for (auto& i0 : closed_) {
            // if this block is not included in the current wave function, skip it

// Some pieces commented out to turn off use of permutation symmetry
//            if (i0.offset() < i2.offset() || i1.offset() < i3.offset()) continue;
            if (!r->is_local(i0, i1, i2, i3) || !r->get_size(i0, i1, i2, i3)) continue;
            const size_t tsize = r->get_size(i0, i1, i2, i3);
            unique_ptr<complex<double>[]> local(new complex<double>[tsize]);
            unique_ptr<complex<double>[]> buf(new complex<double>[tsize]);
            fill_n(buf.get(), tsize, 0.0);

            for (auto& i3t : virt_) {
              for (auto& i1t : virt_) {
//                if (i1t.offset() < i3t.offset()) continue;

                unique_ptr<complex<double>[]> data0 = t->get_block(i0, i1t, i2, i3t);
                unique_ptr<complex<double>[]> data1(new complex<double>[t->get_size(i0, i1t, i2, i3t)]);
                sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1t.size(), i2.size(), i3t.size());

                unique_ptr<complex<double>[]> data2 = v2_->get_block(i1t, i1, i3t, i3);
                unique_ptr<complex<double>[]> data3(new complex<double>[v2_->get_size(i1t, i1, i3t, i3)]);
                sort_indices<0,2,1,3,0,1,1,1>(data2, data3, i1t.size(), i1.size(), i3t.size(), i3.size());

/*
                if (i1t.offset() != i3t.offset()) {
                  unique_ptr<complex<double>[]> data2 = v2_->get_block(i1t, i3, i3t, i1);
                  sort_indices<0,2,3,1,1,1,-1,1>(data2, data3, i1t.size(), i3.size(), i3t.size(), i1.size());
                }
*/

                zgemm3m_("N", "N", i0.size()*i2.size(), i1.size()*i3.size(), i1t.size()*i3t.size(),
                         4.0, data1, i0.size()*i2.size(), data3, i1t.size()*i3t.size(), 1.0, buf, i0.size()*i2.size());
              }
            }
            sort_indices<0,2,1,3,0,1,1,1>(buf, local, i0.size(), i2.size(), i1.size(), i3.size());
            r->add_block(local, i0, i1, i2, i3);
/*
            if (i0.offset() != i2.offset()) {
              sort_indices<2,1,0,3,0,1,-1,1>(local, buf, i0.size(), i1.size(), i2.size(), i3.size());
              r->add_block(buf, i2, i1, i0, i3);
              if (i1.offset() != i3.offset()) {
                sort_indices<0,3,2,1,0,1,-1,1>(buf, local, i2.size(), i1.size(), i0.size(), i3.size());
                r->add_block(local, i2, i3, i0, i1);
                sort_indices<2,1,0,3,0,1,-1,1>(local, buf, i2.size(), i3.size(), i0.size(), i1.size());
                r->add_block(buf, i0, i3, i2, i1);
              }
            } else {
              sort_indices<0,3,2,1,0,1,-1,1>(local, buf, i0.size(), i1.size(), i2.size(), i3.size());
              r->add_block(buf, i0, i3, i2, i1);
            }
*/
          }
        }
      }
    }
  }
  mpi__->barrier();
}

#endif
