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

  public:
    MatView_(const MatView_& o) : btas::TensorView2<DataType>(o), localized_(o.localized()) { }
    MatView_(const btas::TensorView2<DataType>& o, const bool lo) : btas::TensorView2<DataType>(o), localized_(lo) { }
    MatView_(btas::TensorView2<DataType>&& o, const bool lo) : btas::TensorView2<DataType>(std::move(o)), localized_(lo) { }
    MatView_(const btas::CRange<2>& r, const typename btas::Tensor2<DataType>::storage_type& s, const bool lo) : btas::TensorView2<DataType>(r, s), localized_(lo) { }

    int ndim() const { return this->range(0).size(); }
    int mdim() const { return this->range(1).size(); }
    size_t size() const { return ndim() * mdim(); }

    DataType* data()             { assert(contiguous()); return &*this->begin(); }
    const DataType* data() const { assert(contiguous()); return &*this->cbegin(); }

    bool localized() const { return localized_; }
    bool contiguous() const { return this->range().ordinal().contiguous(); }
};

using MatView = MatView_<double>;
using ZMatView = MatView_<std::complex<double>>;

}

#endif
