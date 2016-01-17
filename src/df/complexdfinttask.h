//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexdfinttask.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#ifndef __SRC_DF_COMPLEXDFINTTASK_H
#define __SRC_DF_COMPLEXDFINTTASK_H

#include <src/df/dfblock.h>
#include <src/molecule/shell.h>

namespace bagel {

template <typename TBatch, int N>
class ComplexDFIntTask {
  protected:
    const std::array<std::shared_ptr<const Shell>,4> shell_;
    const std::array<int,3> offset_; // at most 3 elements
    std::array<std::shared_ptr<DFBlock>,N> dfblocks_;

    std::shared_ptr<TBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,4>& input) const {
      auto eribatch = std::make_shared<TBatch>(input, 2.0);
      assert(eribatch->Nblocks() * 2 == N); // Twice as many blocks so we can separate real and imaginary parts
      eribatch->compute();
      return eribatch;
    }

  public:
    ComplexDFIntTask(std::array<std::shared_ptr<const Shell>,4>&& a, std::array<int,3>&& b, std::array<std::shared_ptr<DFBlock>,N>& df)
     : shell_(a), offset_(b), dfblocks_(df) { };

    void compute() {
      std::shared_ptr<TBatch> p = compute_batch(shell_);
      const int Nhalf = N / 2;

      // all slot in
      for (int i = 0; i != Nhalf; ++i) {
        const int i2 = 2*i;
        assert(dfblocks_[i2]->b1size() == dfblocks_[i2]->b2size());
        assert(dfblocks_[i2]->asize() == dfblocks_[i2+1]->asize());
        assert(dfblocks_[i2]->b1size() == dfblocks_[i2+1]->b1size());
        assert(dfblocks_[i2]->b2size() == dfblocks_[i2+1]->b2size());
        const size_t nbin = dfblocks_[i2]->b1size();
        const size_t naux = dfblocks_[i2]->asize();

        const std::complex<double>* ppt = p->data(i);
        double* const data_r = dfblocks_[i2]->data();
        double* const data_i = dfblocks_[i2+1]->data();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1, ppt += shell_[1]->nbasis()) {
            for (int n=0; n!=shell_[1]->nbasis(); ++n) data_r[offset_[2]+naux*(j1+nbin*j0)+n] = std::real(ppt[n]);
            for (int n=0; n!=shell_[1]->nbasis(); ++n) data_i[offset_[2]+naux*(j1+nbin*j0)+n] = std::imag(ppt[n]);
            if (N == 2) {
              for (int n=0; n!=shell_[1]->nbasis(); ++n) data_r[offset_[2]+naux*(j0+nbin*j1)+n] = std::real(ppt[n]);
              for (int n=0; n!=shell_[1]->nbasis(); ++n) data_i[offset_[2]+naux*(j0+nbin*j1)+n] = -std::imag(ppt[n]);
            }
          }
        }
      }
    }
};


}

#endif
