//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rmawindow.h
// Copyright (C) 2016 Toru Shiozaki
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

#ifndef __SRC_PARALLEL_RMAWINDOW_H
#define __SRC_PARALLEL_RMAWINDOW_H

#include <bagel_config.h>
#include <complex>
#include <memory>
#ifdef HAVE_MPI_H
 #include <mpi.h>
#endif

namespace bagel {

// task class
template<typename DataType>
class RMATask {
  public:
#ifndef HAVE_MPI_H
    using MPI_Request = int; // just to compile
#endif
  protected:
    MPI_Request tag;
    std::unique_ptr<DataType[]> buf;
  public:
    RMATask(MPI_Request&& t) : tag(std::move(t)) { }
    RMATask(MPI_Request&& t, std::unique_ptr<DataType[]>&& o) : tag(std::move(t)), buf(std::move(o)) { }
    void wait();
    bool test();
};


// window class
template<typename DataType>
class RMAWindow {
  protected:
#ifndef HAVE_MPI_H
    using MPI_Win = int; // just to compile
#endif
    MPI_Win win_;
    DataType* win_base_;

    bool initialized_;

  public:
    RMAWindow();
    virtual ~RMAWindow();

    RMAWindow<DataType>& operator=(const RMAWindow<DataType>& o);

    void initialize();
    bool initialized() const { return initialized_; }
    void zero();
    void scale(const DataType& a);

    void fence() const;
    void fence_local() const;

    void ax_plus_y(const DataType& a, const RMAWindow<DataType>& o);
    void ax_plus_y(const DataType& a, std::shared_ptr<const RMAWindow<DataType>> o) { ax_plus_y(a, *o); }

    DataType dot_product(const RMAWindow<DataType>& o) const;
    DataType dot_product(std::shared_ptr<const RMAWindow<DataType>> o) const { return dot_product(*o); }

    const DataType* local_data() const { fence(); return win_base_; }

    // Blocking
    std::unique_ptr<DataType[]> rma_get(const size_t key) const;
    std::unique_ptr<DataType[]> rma_get(const size_t rank, const size_t off, const size_t size) const;
    void rma_get(DataType*, const size_t key) const;
    void rma_get(DataType*, const size_t rank, const size_t off, const size_t size) const;
    void rma_put(const std::unique_ptr<DataType[]>& dat, const size_t key) { rma_put(dat.get(), key); }
    void rma_put(const std::unique_ptr<DataType[]>& dat, const size_t rank, const size_t off, const size_t size) { rma_put(dat.get(), rank, off, size); }
    void rma_put(const DataType* dat, const size_t key);
    void rma_put(const DataType* dat, const size_t rank, const size_t off, const size_t size);
    void rma_add(const std::unique_ptr<DataType[]>& dat, const size_t key);
    void rma_add(const std::unique_ptr<DataType[]>& dat, const size_t rank, const size_t off, const size_t size) { rma_add(dat.get(), rank, off, size); }
    void rma_add(const DataType* dat, const size_t rank, const size_t off, const size_t size);

    // Non-blocking: The buffer is pushed to the RMATask
    std::shared_ptr<RMATask<DataType>> rma_radd(std::unique_ptr<DataType[]>&& dat, const size_t key);
    std::shared_ptr<RMATask<DataType>> rma_radd(std::unique_ptr<DataType[]>&& dat, const size_t rank, const size_t off, const size_t size);
    // Non-blocking: User needs to manage memory
    std::shared_ptr<RMATask<DataType>> rma_rget(DataType* dat, const size_t key) const;
    std::shared_ptr<RMATask<DataType>> rma_rget(DataType* dat, const size_t rank, const size_t off, const size_t size) const;
    std::shared_ptr<RMATask<DataType>> rma_rput(const DataType* dat, const size_t rank, const size_t off, const size_t size);
    std::shared_ptr<RMATask<DataType>> rma_radd(const DataType* dat, const size_t rank, const size_t off, const size_t size);

    void set_element(const size_t rank, const size_t disp, const DataType a);

    // this has to be called collectively
    void accumulate_buffer(const DataType a, const std::unique_ptr<DataType[]>& buf);

    // returns (process, offset, size)
    virtual bool is_local(const size_t key) const = 0;
    virtual std::tuple<size_t, size_t, size_t> locate(const size_t key) const = 0;
    virtual size_t localsize() const = 0;
};


template<typename DataType>
class RMAWindow_bare : public RMAWindow<DataType> {
  protected:
    size_t localsize_;
  public:
    RMAWindow_bare(const size_t size) : localsize_(size) { this->initialize(); }

    bool is_local(const size_t key) const { throw std::logic_error("RMAWindow_bare::is_local should not be called"); }
    std::tuple<size_t, size_t, size_t> locate(const size_t key) const { throw std::logic_error("RMAWindow_bare::locate should not be called"); }
    size_t localsize() const override { return localsize_; }
};


extern template class RMATask<double>;
extern template class RMATask<std::complex<double>>;
extern template class RMAWindow<double>;
extern template class RMAWindow<std::complex<double>>;

}

#endif
