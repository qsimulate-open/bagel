//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: preallocarray.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_MATH_PREALLOCARRAY_H
#define __SRC_MATH_PREALLOCARRAY_H

#include <complex>
#include <stdexcept>

namespace bagel {

template <class DataType>
class PreAllocArray {
  public:
    using value_type = DataType;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using iterator = DataType*;
    using reference = DataType&;

    using const_value_type = const typename std::remove_cv<DataType>::type;
    using const_iterator = const typename std::remove_cv<DataType>::type*;
    using const_reference = const typename std::remove_cv<DataType>::type&;

  protected:
    DataType* begin_;
    DataType* end_;

  public:
    PreAllocArray() : begin_(0), end_(0) { }
    PreAllocArray(DataType* b, size_t size) : begin_(b), end_(b+size) { }
    PreAllocArray(DataType* b, DataType* e) : begin_(b), end_(e) { }
    PreAllocArray(const PreAllocArray<DataType>& o) : begin_(o.begin_), end_(o.end_) { }

    template <typename T = DataType, class = typename std::enable_if<not std::is_const<T>::value>::type>
    iterator begin() { return begin_; }
    const_iterator begin() const { return begin_; }
    const_iterator cbegin() const { return begin_; }

    template <typename T = DataType, class = typename std::enable_if<not std::is_const<T>::value>::type>
    iterator end() { return end_; }
    const_iterator end() const { return end_; }
    const_iterator cend() const { return end_; }

    size_t size() const { return std::distance(begin_, end_); }

    void resize(const int) { throw std::logic_error("resize is not allowed in PreAllocArray"); }
    bool empty() const { return begin_ == end_; }

    template <typename T = DataType, class = typename std::enable_if<not std::is_const<T>::value>::type>
    value_type& at(const size_t i) { return *(begin_+i); }
    const_value_type& at(const size_t i) const { return *(begin_+i); }

    template <typename T = DataType, class = typename std::enable_if<not std::is_const<T>::value>::type>
    value_type& front() { return at(0); }
    const_value_type& front() const { return at(0); }
    template <typename T = DataType, class = typename std::enable_if<not std::is_const<T>::value>::type>
    value_type& back() { return at(size()); }
    const_value_type& back() const { return at(size()); }

    template <typename T = DataType, class = typename std::enable_if<not std::is_const<T>::value>::type>
    iterator data() { return begin_; }
    const_iterator data() const { return begin_; }

    PreAllocArray<DataType>& operator=(const PreAllocArray<DataType>& o) { begin_ = o.begin_; end_ = o.end_; return *this; }
};

}

extern template class bagel::PreAllocArray<double>;
extern template class bagel::PreAllocArray<std::complex<double>>;
extern template class bagel::PreAllocArray<const double>;
extern template class bagel::PreAllocArray<const std::complex<double>>;

#endif
