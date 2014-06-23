//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock_base.h
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

#ifndef __SRC_DF_DFBLOCK_BASE_H
#define __SRC_DF_DFBLOCK_BASE_H

#include <numeric>
#include <src/math/algo.h>
#include <src/util/timer.h>
#include <src/util/simple.h>
#include <src/util/taskqueue.h>
#include <src/parallel/staticdist.h>
#include <src/parallel/mpi_interface.h>
#include <src/math/btas_interface.h>

namespace bagel {

template <typename DataType>
class DFBlock_base : public btas::Tensor3<DataType> {
  // aux_ runs fastest, b2_ runs slowest
  public:
    using btas::Tensor3<DataType>::data;

  protected:
    // distribution information
    std::shared_ptr<const StaticDist> adist_shell_;
    std::shared_ptr<const StaticDist> adist_;
    // if true, asize is evenly distributed. If false, asize is at the shell boundary
    bool averaged_;

    // a set of offsets of this block in the entire DF integrals
    size_t astart_;
    size_t b1start_;
    size_t b2start_;

  public:
    DFBlock_base() { }

    // construction of a block from AO integrals
    DFBlock_base(std::shared_ptr<const StaticDist> adist_shell, std::shared_ptr<const StaticDist> adist,
                 const size_t a, const size_t b1, const size_t b2, const int as, const int b1s, const int b2s, const bool averaged = false)
     : btas::Tensor3<DataType>(std::max(adist_shell->size(mpi__->rank()), std::max(adist->size(mpi__->rank()), a)), b1, b2),
       adist_shell_(adist_shell), adist_(adist), averaged_(averaged), astart_(as), b1start_(b1s), b2start_(b2s) {

      assert(asize() == adist_shell->size(mpi__->rank()) || asize() == adist_->size(mpi__->rank()) || asize() == adist_->nele());

      // resize to the current size (moving the end pointer)
      const btas::CRange<3> range(a, b1, b2);
      this->resize(range);
    }

    DFBlock_base(const DFBlock_base<DataType>& o)
     : btas::Tensor3<DataType>(std::max(o.adist_shell_->size(mpi__->rank()), std::max(o.adist_->size(mpi__->rank()), o.asize())), o.b1size(), o.b2size()),
       adist_shell_(o.adist_shell_), adist_(o.adist_), averaged_(o.averaged_), astart_(o.astart_), b1start_(o.b1start_), b2start_(o.b2start_) {

      // resize to the current size
      const btas::CRange<3> range(o.asize(), o.b1size(), o.b2size());
      this->resize(range);

      btas::Tensor3<DataType>::operator=(o);
    }

    // dimensions of the block
    size_t asize() const { return this->extent(0); }
    size_t b1size() const { return this->extent(1); }
    size_t b2size() const { return this->extent(2); }

    size_t size() const { return asize()*b1size()*b2size(); }
    bool averaged() const { return averaged_; }

    // a set of offsets of this block in the entire DF integrals
    size_t astart() const { return astart_; }
    size_t b1start() const { return b1start_; }
    size_t b2start() const { return b2start_; }

    // dist
    const std::shared_ptr<const StaticDist>& adist_now() const { return averaged_ ? adist_ : adist_shell_; }

    // some math functions
    DFBlock_base<DataType>& operator=(const DFBlock_base<DataType>& o) {
      btas::Tensor3<DataType>::operator=(o);
      adist_shell_ = o.adist_shell_;
      adist_ = o.adist_;
      averaged_ = o.averaged_;
      astart_ = o.astart_;
      b1start_ = o.b1start_;
      b2start_ = o.b2start_;
      return *this;
    }
    DFBlock_base<DataType>& operator=(DFBlock_base<DataType>&& o) {
      btas::Tensor3<DataType>::operator=(std::move(o));
      adist_shell_ = o.adist_shell_;
      adist_ = o.adist_;
      averaged_ = o.averaged_;
      astart_ = o.astart_;
      b1start_ = o.b1start_;
      b2start_ = o.b2start_;
      return *this;
    }
    DFBlock_base<DataType>& operator+=(const DFBlock_base<DataType>& o) { btas::Tensor3<DataType>::operator+=(o); return *this; }
    DFBlock_base<DataType>& operator-=(const DFBlock_base<DataType>& o) { btas::Tensor3<DataType>::operator-=(o); return *this; }

    template <typename ScaleType, class DType>
    void ax_plus_y(const ScaleType a, const DType& o) { btas::axpy(a, o, *this); }
    template <typename ScaleType, class DType>
    void ax_plus_y(const ScaleType a, const std::shared_ptr<DType>& o) { ax_plus_y(a, *o); }
    template <typename ScaleType>
    void scale(const ScaleType a) { btas::scal(a, *this); }

    void zero() { std::fill_n(data(), size(), 0.0); }

    // symmetrize b1 and b2 (assuming b1size() == b2size())
    void symmetrize() {
      if (b1size() != b2size()) throw std::logic_error("illegal call of DFBlock::symmetrize()");
      const int n = b1size();
      for (int i = 0; i != n; ++i)
        for (int j = i; j != n; ++j) {
          blas::ax_plus_y_n(1.0, data()+asize()*(j+n*i), asize(), data()+asize()*(i+n*j));
          std::copy_n(data()+asize()*(i+n*j), asize(), data()+asize()*(j+n*i));
        }
    }

    template <class MatType>
    void copy_block(std::shared_ptr<MatType> o, const int jdim, const size_t offset) {
      assert(o->size() == asize()*jdim);
      std::copy_n(o->data(), asize()*jdim, data()+offset);
    }
    template <class MatType>
    void add_block(std::shared_ptr<MatType> o, const int jdim, const size_t offset, const DataType fac = 1.0) {
      assert(o->size() == asize()*jdim);
      blas::ax_plus_y_n(fac, o->data(), asize()*jdim, data()+offset);
    }


    // average the asize between MPI processes (block will be described by dist_)
    void average() {
      if (averaged_) return;
      averaged_ = true;

      // first make a send and receive buffer
      const size_t o_start = astart_;
      const size_t o_end   = o_start + asize();
      const int myrank = mpi__->rank();
      size_t t_start, t_end;
      std::tie(t_start, t_end) = adist_->range(myrank);

      assert(o_end - t_end >= 0);
      assert(o_start - t_start >= 0);

      // TODO so far I am not considering the cases when data must be sent to the next neighbor; CAUTION
      const size_t asendsize = o_end - t_end;
      const size_t arecvsize = o_start - t_start;

      assert(asendsize < t_end-t_start && arecvsize < t_end-t_start);

      std::unique_ptr<DataType[]> sendbuf;
      std::unique_ptr<DataType[]> recvbuf;
      int sendtag = 0;
      int recvtag = 0;

      if (asendsize) {
        TaskQueue<CopyBlockTask<DataType>> task(b2size());

        sendbuf = std::unique_ptr<DataType[]>(new DataType[asendsize*b1size()*b2size()]);
        const size_t retsize = asize() - asendsize;
        for (size_t b2 = 0; b2 != b2size(); ++b2)
          task.emplace_back(data()+retsize+asize()*b1size()*b2, asize(), sendbuf.get()+asendsize*b1size()*b2, asendsize, asendsize, b1size());

        task.compute();

        // send to the next node
        sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size()*b2size(), myrank+1, myrank);
      }

      if (arecvsize) {
        recvbuf = std::unique_ptr<DataType[]>(new DataType[arecvsize*b1size()*b2size()]);
        // recv from the previous node
        recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size()*b2size(), myrank-1, myrank-1);
      }

      // second move local data
      if (arecvsize || asendsize) {
        const size_t t_size = t_end - t_start;
        const size_t retsize = asize() - asendsize;
        if (t_size <= asize()) {
          for (size_t i = 0; i != b1size()*b2size(); ++i) {
            if (i*asize() < (i+1)*t_size-retsize) {
              std::copy_backward(data()+i*asize(), data()+i*asize()+retsize, data()+(i+1)*t_size);
            } else if (i*asize() > (i+1)*t_size-retsize) {
              std::copy_n(data()+i*asize(), retsize, data()+(i+1)*t_size-retsize);
            }
          }
        } else {
          for (long long int i = b1size()*b2size()-1; i >= 0; --i) {
            assert(i*asize() < (i+1)*t_size-retsize);
            std::copy_backward(data()+i*asize(), data()+i*asize()+retsize, data()+(i+1)*t_size);
          }
        }
      }

      // set new astart_ and asize()
      astart_ = t_start;
      assert(this->storage().capacity() >= (t_end - t_start)*b1size()*b2size());
      const btas::CRange<3> range(t_end - t_start, b1size(), b2size());
      this->resize(range);

      // set received data
      if (arecvsize) {
        // wait for recv communication
        mpi__->wait(recvtag);

        TaskQueue<CopyBlockTask<DataType>> task(b2size());
        for (size_t b2 = 0; b2 != b2size(); ++b2)
          task.emplace_back(recvbuf.get()+arecvsize*b1size()*b2, arecvsize, data()+asize()*b1size()*b2, asize(), arecvsize, b1size());
        task.compute();
      }

      // wait for send communication
      if (asendsize) mpi__->wait(sendtag);
    }

    // reverse operation of average() function
    void shell_boundary() {
      if (!averaged_) return;
      averaged_ = false;
      const size_t o_start = astart_;
      const size_t o_end = o_start + asize();
      const int myrank = mpi__->rank();
      size_t t_start, t_end;
      std::tie(t_start, t_end) = adist_shell_->range(myrank);

      const size_t asendsize = t_start - o_start;
      const size_t arecvsize = t_end - o_end;
      assert(t_start >= o_start && t_end >= o_end);

      std::unique_ptr<DataType[]> sendbuf, recvbuf;
      int sendtag = 0;
      int recvtag = 0;

      if (asendsize) {
        TaskQueue<CopyBlockTask<DataType>> task(b2size());
        sendbuf = std::unique_ptr<DataType[]>(new DataType[asendsize*b1size()*b2size()]);
        for (size_t b2 = 0; b2 != b2size(); ++b2)
          task.emplace_back(data()+asize()*b1size()*b2, asize(), sendbuf.get()+asendsize*b1size()*b2, asendsize, asendsize, b1size());

        task.compute();
        assert(myrank > 0);
        sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size()*b2size(), myrank-1, myrank);
      }
      if (arecvsize) {
        assert(myrank+1 < mpi__->size());
        recvbuf = std::unique_ptr<DataType[]>(new DataType[arecvsize*b1size()*b2size()]);
        recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size()*b2size(), myrank+1, myrank+1);
      }

      if (arecvsize || asendsize) {
        const size_t t_size = t_end - t_start;
        const size_t retsize = asize() - asendsize;
        assert(t_size >= retsize);
        if (t_size <= asize()) {
          for (size_t i = 0; i != b1size()*b2size(); ++i) {
            assert(i*asize()+asendsize > i*t_size);
            std::copy_n(data()+i*asize()+asendsize, retsize, data()+i*t_size);
          }
        } else {
          for (long long int i = b1size()*b2size()-1; i >= 0; --i) {
            if (i*asize()+asendsize > i*t_size) {
              std::copy_n(data()+i*asize()+asendsize, retsize, data()+i*t_size);
            } else if (i*asize()+asendsize < i*t_size) {
              std::copy_backward(data()+i*asize()+asendsize, data()+(i+1)*asize(), data()+i*t_size+retsize);
            }
          }
        }
      }

      // set new astart_ and asize()
      astart_ = t_start;
      assert(this->storage().capacity() >= (t_end - t_start)*b1size()*b2size());
      const btas::CRange<3> range(t_end - t_start, b1size(), b2size());
      this->resize(range);

      // set received data
      if (arecvsize) {
        // wait for recv communication
        mpi__->wait(recvtag);

        TaskQueue<CopyBlockTask<DataType>> task(b2size());
        for (size_t b2 = 0; b2 != b2size(); ++b2)
          task.emplace_back(recvbuf.get()+arecvsize*b1size()*b2, arecvsize, data()+asize()*b1size()*b2+(asize()-arecvsize), asize(), arecvsize, b1size());
        task.compute();
      }

      // wait for send communication
      if (asendsize) mpi__->wait(sendtag);
    }
};

extern template class DFBlock_base<double>;
extern template class DFBlock_base<std::complex<double>>;

}

#endif
