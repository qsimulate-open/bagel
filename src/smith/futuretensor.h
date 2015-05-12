//
// BAGEL - Parallel electron correlation program.
// Filename: futuretensor.h
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


#ifndef __SRC_SMITH_FUTURETENSOR_H
#define __SRC_SMITH_FUTURETENSOR_H

#include <src/smith/tensor.h>
#include <src/smith/task.h>

namespace bagel {
namespace SMITH {

class FutureTensor : public Tensor {
  protected:
    // TODO actually not const, but this is the only way to make it compiled...
    void init() const override { init_->compute(); initialized_ = true; }
    mutable std::shared_ptr<Task> init_;

  public:
    FutureTensor(const Tensor& i,  std::shared_ptr<Task> j) : Tensor(i), init_(j) { }

};

}
}

#endif
