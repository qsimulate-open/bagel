//
// BAGEL - Parallel electron correlation program.
// Filename: preallocarray.h
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

#ifndef __SRC_MATH_PREALLOCARRAY_H
#define __SRC_MATH_PREALLOCARRAY_H

#include <complex>
#include <stdexcept>

namespace bagel {

template <class DataType>
class PreAllocArray_ {
  public:
    using value_type = DataType;
    using size_type = size_t;
    using iterator = DataType*;
    using const_iterator = const DataType*;

  protected:
    DataType* begin_;
    DataType* end_;

  public:
    PreAllocArray_() : begin_(0), end_(0) { }
    PreAllocArray_(DataType* b, size_t size) : begin_(b), end_(b+size) { }
    PreAllocArray_(DataType* b, DataType* e) : begin_(b), end_(e) { }
    PreAllocArray_(const PreAllocArray_<DataType>& o) : begin_(o.begin_), end_(o.end_) { }

    DataType* begin() { return begin_; }
    const DataType* begin() const { return begin_; }
    const DataType* cbegin() const { return begin_; }

    DataType* end() { return end_; }
    const DataType* end() const { return end_; }
    const DataType* cend() const { return end_; }

    size_t size() const { return std::distance(begin_, end_); }

    void resize(const int) { throw std::logic_error("resize is not allowed in PreAllocArray_"); }
    bool empty() const { return begin_ == end_; }

    DataType& at(const size_t i) { return *(begin_+i); }
    const DataType& at(const size_t i) const { return *(begin_+i); }

    DataType& front() { return at(0); }
    const DataType& front() const { return at(0); }
    DataType& back() { return at(size()); }
    const DataType& back() const { return at(size()); }

    DataType* data() { return begin_; }
    const DataType* data() const { return begin_; }

    PreAllocArray_<DataType>& operator=(const PreAllocArray_<DataType>& o) { begin_ = o.begin_; end_ = o.end_; return *this; }
};

using PreAllocArray = PreAllocArray_<double>;
using ZPreAllocArray = PreAllocArray_<std::complex<double>>;

}

extern template class bagel::PreAllocArray_<double>;
extern template class bagel::PreAllocArray_<std::complex<double>>;

#endif
