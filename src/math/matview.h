//
// BAGEL - Parallel electron correlation program.
// Filename: matview.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_MATH_MATVIEW_H
#define __SRC_MATH_MATVIEW_H

#include <complex>
#include <src/math/btas_interface.h>

namespace bagel {

template<typename DataType>
class MatView_ : public btas::TensorView2<DataType> {
  protected:
    bool localized_;

#ifdef HAVE_SCALAPACK
    std::vector<int> desc_;
    std::tuple<int, int> localsize_;

    void setlocal_(const std::unique_ptr<DataType[]>& local) {
      zero();

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
            std::copy_n(&local[j*blocksize__+localrow*(i*blocksize__+id)], blocksize__, element_ptr(myprow+j*nstride, mypcol+i*mstride+id));

      for (int id = 0; id != localcol % blocksize__; ++id) {
        for (int j = 0; j != nblock; ++j)
          std::copy_n(&local[j*blocksize__+localrow*(mblock*blocksize__+id)], blocksize__, element_ptr(myprow+j*nstride, mypcol+mblock*mstride+id));
        for (int jd = 0; jd != localrow % blocksize__; ++jd)
          element(myprow+nblock*nstride+jd, mypcol+mblock*mstride+id) = local[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)];
      }
      for (int i = 0; i != mblock; ++i)
        for (int id = 0; id != blocksize__; ++id)
          for (int jd = 0; jd != localrow % blocksize__; ++jd)
            element(myprow+nblock*nstride+jd, mypcol+i*mstride+id) = local[nblock*blocksize__+jd+localrow*(i*blocksize__+id)];

      // syncronize (this can be improved, but...)
      allreduce();
    }
#endif

    void init() {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim(), mdim());
        localsize_ = mpi__->numroc(ndim(), mdim());
      }
#endif
    }

  public:
    MatView_(const MatView_& o) : btas::TensorView2<DataType>(o), localized_(o.localized()) { init(); }
    MatView_(const btas::TensorView2<DataType>& o, const bool lo) : btas::TensorView2<DataType>(o), localized_(lo) { init(); }
    MatView_(btas::TensorView2<DataType>&& o, const bool lo) : btas::TensorView2<DataType>(std::move(o)), localized_(lo) { init(); }
    MatView_(const btas::CRange<2>& r, const typename btas::Tensor2<DataType>::storage_type& s, const bool lo) : btas::TensorView2<DataType>(r, s), localized_(lo) { init(); }

    int ndim() const { return this->extent(0); }
    int mdim() const { return this->extent(1); }
    size_t size() const { return ndim() * mdim(); }

    DataType* data()             { assert(contiguous()); return &*this->begin(); }
    const DataType* data() const { assert(contiguous()); return &*this->cbegin(); }

    bool localized() const { return localized_; }
    bool contiguous() const { return this->range().ordinal().contiguous(); }

    void zero() const { std::fill_n(this->begin(), this->end(), 0.0); }
    DataType& element(const int i, const int j) { return *element_ptr(i,j); }
    const DataType& element(const int i, const int j) const { return *element_ptr(i,j); }
    DataType* element_ptr(const int i, const int j) { return data()+i+ndim()*j; }
    const DataType* element_ptr(const int i, const int j) const { return data()+i+ndim()*j; }

    void allreduce() {
      assert(!localized_);
      mpi__->allreduce(data(), size());
    }

#ifdef HAVE_SCALAPACK
    const std::vector<int>& desc() const { return desc_; }
    void setlocal(const std::unique_ptr<DataType[]>& local) { setlocal_(local); }

    std::unique_ptr<DataType[]> getlocal() const {
      const int localrow = std::get<0>(localsize_);
      const int localcol = std::get<1>(localsize_);

      std::unique_ptr<DataType[]> local(new DataType[localrow*localcol]);

      const int nblock = localrow/blocksize__;
      const int mblock = localcol/blocksize__;
      const size_t nstride = blocksize__*mpi__->nprow();
      const size_t mstride = blocksize__*mpi__->npcol();
      const int myprow = mpi__->myprow()*blocksize__;
      const int mypcol = mpi__->mypcol()*blocksize__;

      for (int i = 0; i != mblock; ++i)
        for (int j = 0; j != nblock; ++j)
          for (int id = 0; id != blocksize__; ++id)
            std::copy_n(element_ptr(myprow+j*nstride, mypcol+i*mstride+id), blocksize__, &local[j*blocksize__+localrow*(i*blocksize__+id)]);

      for (int id = 0; id != localcol % blocksize__; ++id) {
        for (int j = 0; j != nblock; ++j)
          std::copy_n(element_ptr(myprow+j*nstride, mypcol+mblock*mstride+id), blocksize__, &local[j*blocksize__+localrow*(mblock*blocksize__+id)]);
        for (int jd = 0; jd != localrow % blocksize__; ++jd)
          local[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)] = element(myprow+nblock*nstride+jd, mypcol+mblock*mstride+id);
      }
      for (int i = 0; i != mblock; ++i)
        for (int id = 0; id != blocksize__; ++id)
          for (int jd = 0; jd != localrow % blocksize__; ++jd)
            local[nblock*blocksize__+jd+localrow*(i*blocksize__+id)] = element(myprow+nblock*nstride+jd, mypcol+i*mstride+id);
      return local;
    }
#endif
};

using MatView = MatView_<double>;
using ZMatView = MatView_<std::complex<double>>;

}

#endif
