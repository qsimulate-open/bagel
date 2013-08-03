//
// BAGEL - Parallel electron correlation program.
// Filename: distmatrix_base.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_UTIL_DISTMATRIX_BASE_H
#define __SRC_UTIL_DISTMATRIX_BASE_H

#include <cassert>
#include <algorithm>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

namespace bagel {

#ifdef HAVE_SCALAPACK

template<typename DataType>
class DistMatrix_base {
  protected:
    // global dimension
    const int ndim_;
    const int mdim_;

    // distributed data
    std::unique_ptr<DataType[]> local_;

    // Scalapack specific
    const std::unique_ptr<int[]> desc_;
    const std::tuple<int, int> localsize_;

  public:
    DistMatrix_base(const int n, const int m) : ndim_(n), mdim_(m), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
      local_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      zero();
    }
    DistMatrix_base(const DistMatrix_base& o) : ndim_(o.ndim_), mdim_(o.mdim_), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
      local_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      std::copy_n(o.local_.get(), size(), local_.get());
    }

    const std::unique_ptr<DataType[]>& local() const { return local_; }

    size_t size() const { return std::get<0>(localsize_)*std::get<1>(localsize_); }
    int ndim() const { return ndim_; }
    int mdim() const { return mdim_; }

    virtual void diagonalize(double* vec) = 0;

    void fill(const DataType a) { std::fill_n(local_.get(), size(), a); }
    void zero() { const DataType zero(0.0); fill(zero); }

    void add_diag(const DataType& a, const size_t start, const size_t fence) {
      assert(ndim_ == mdim_ && start <= fence);
      const int localrow = std::get<0>(localsize_);
      const int localcol = std::get<1>(localsize_);

      const int nblock = localrow/blocksize__;
      const int mblock = localcol/blocksize__;
      const size_t nstride = blocksize__*mpi__->nprow();
      const size_t mstride = blocksize__*mpi__->npcol();
      const int myprow = mpi__->myprow()*blocksize__;
      const int mypcol = mpi__->mypcol()*blocksize__;

      for (int i = 0; i != mblock; ++i)
        for (int j = 0; j != nblock; ++j)
          for (int id = 0; id != blocksize__; ++id)
            for (int jd = 0; jd != blocksize__; ++jd) {
              const size_t m = mypcol+i*mstride+id;
              const size_t n = myprow+j*nstride+jd;
              if (m == n && start <= n && n < fence)
                local_[j*blocksize__+jd+localrow*(i*blocksize__+id)] += a;
            }
        for (int id = 0; id != localcol % blocksize__; ++id)
          for (int jd = 0; jd != localrow % blocksize__; ++jd) {
            const size_t m = mblock*blocksize__+id;
            const size_t n = nblock*blocksize__+jd;
            if (m == n && start <= n && n < fence)
              local_[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)] += a;
          }
    }

    void scale(const double* vec) {
      const int localrow = std::get<0>(localsize_);
      const int localcol = std::get<1>(localsize_);

      const int nblock = localrow/blocksize__;
      const int mblock = localcol/blocksize__;
      const size_t nstride = blocksize__*mpi__->nprow();
      const size_t mstride = blocksize__*mpi__->npcol();
      const int myprow = mpi__->myprow()*blocksize__;
      const int mypcol = mpi__->mypcol()*blocksize__;
      for (int i = 0; i != mblock; ++i)
        for (int j = 0; j != nblock; ++j)
          for (int id = 0; id != blocksize__; ++id) {
            DataType* c = local_.get()+j*blocksize__+localrow*(i*blocksize__+id);
            std::for_each(c, c+blocksize__, [&x = vec[mypcol+i*mstride+id]](DataType& a) { a*= x; });
          }

      for (int id = 0; id != localcol % blocksize__; ++id) {
        for (int j = 0; j != nblock; ++j) {
          DataType* c = local_.get()+j*blocksize__+localrow*(mblock*blocksize__+id);
          std::for_each(c, c+blocksize__, [&x = vec[mypcol+mblock*mstride+id]](DataType& a) { a*= x; });
        }

        DataType* c = local_.get()+nblock*blocksize__+localrow*(mblock*blocksize__+id);
        std::for_each(c, c+(localrow%blocksize__), [&x = vec[mypcol+mblock*mstride+id]](DataType& a) { a*= x; });
      }
      for (int i = 0; i != mblock; ++i)
        for (int id = 0; id != blocksize__; ++id) {
          DataType* c = local_.get()+nblock*blocksize__+localrow*(i*blocksize__+id);
          std::for_each(c, c+(localrow%blocksize__), [&x = vec[mypcol+i*mstride+id]](DataType& a) { a*= x; });
        }
    }
};
#endif

}

#endif
