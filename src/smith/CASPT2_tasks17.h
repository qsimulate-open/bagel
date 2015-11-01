//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks17.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#ifndef __SRC_SMITH_CASPT2_TASKS17_H
#define __SRC_SMITH_CASPT2_TASKS17_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task800 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task800(std::shared_ptr<TATensor<double,2>> I1082, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1082), ta1_(t2) { }
};

class Task801 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task801(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma362, std::shared_ptr<TATensor<double,4>> I1124)
   : ta0_(I768), ta1_(Gamma362), ta2_(I1124) { }
};

class Task802 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task802(std::shared_ptr<TATensor<double,4>> I1124, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1124), ta1_(t2) { }
};

class Task803 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task803(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma377, std::shared_ptr<TATensor<double,6>> I1170)
   : ta0_(I768), ta1_(Gamma377), ta2_(I1170) { }
};

class Task804 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task804(std::shared_ptr<TATensor<double,6>> I1170, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1170), ta1_(v2), ta2_(t2) { }
};

class Task805 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task805(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma395, std::shared_ptr<TATensor<double,6>> I1224)
   : ta0_(I768), ta1_(Gamma395), ta2_(I1224) { }
};

class Task806 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task806(std::shared_ptr<TATensor<double,6>> I1224, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1224), ta1_(v2), ta2_(t2) { }
};


}
}
}
#endif
#endif

