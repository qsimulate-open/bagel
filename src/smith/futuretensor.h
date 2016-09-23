//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: futuretensor.h
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


#ifndef __SRC_SMITH_FUTURETENSOR_H
#define __SRC_SMITH_FUTURETENSOR_H

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/tensor.h>
#include <src/smith/task.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class FutureTensor_ : public Tensor_<DataType> {
  protected:
    using Tensor_<DataType>::initialized_;

  protected:
    mutable std::shared_ptr<Task> init_;

  public:
    FutureTensor_(const Tensor_<DataType>& i,  std::shared_ptr<Task> j);

    void init() const override;

};

extern template class FutureTensor_<double>;
extern template class FutureTensor_<std::complex<double>>;

namespace CASPT2 { using FutureTensor = FutureTensor_<double>; }
namespace CASA { using FutureTensor = FutureTensor_<double>; }
namespace SPCASPT2 { using FutureTensor = FutureTensor_<double>; }
namespace MSCASPT2 { using FutureTensor = FutureTensor_<double>; }
namespace MRCI   { using FutureTensor = FutureTensor_<double>; }
namespace RelCASPT2 { using FutureTensor = FutureTensor_<std::complex<double>>; }
namespace RelCASA { using FutureTensor = FutureTensor_<std::complex<double>>; }
namespace RelMRCI   { using FutureTensor = FutureTensor_<std::complex<double>>; }

}
}

#endif
#endif
