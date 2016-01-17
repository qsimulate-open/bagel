//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks11.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#ifndef __SRC_SMITH_CASPT2_TASKS11_H
#define __SRC_SMITH_CASPT2_TASKS11_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task500 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task500(std::shared_ptr<TATensor<double,4>> I644, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I644), ta1_(Gamma38), ta2_(t2) { }
};

class Task501 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task501(std::shared_ptr<TATensor<double,2>> I643, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I647)
   : ta0_(I643), ta1_(t2), ta2_(I647) { }
};

class Task502 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task502(std::shared_ptr<TATensor<double,4>> I647, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I647), ta1_(Gamma38), ta2_(t2) { }
};

class Task503 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task503(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I649)
   : ta0_(proj), ta1_(I649) { }
};

class Task504 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task504(std::shared_ptr<TATensor<double,2>> I649, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I650)
   : ta0_(I649), ta1_(t2), ta2_(I650) { }
};

class Task505 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task505(std::shared_ptr<TATensor<double,4>> I650, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I650), ta1_(Gamma60), ta2_(t2) { }
};

class Task506 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task506(std::shared_ptr<TATensor<double,2>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task507 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task507(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I664)
   : ta0_(proj), ta1_(I664) { }
};

class Task508 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task508(std::shared_ptr<TATensor<double,2>> I664, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I664), ta1_(Gamma12), ta2_(t2) { }
};

class Task509 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task509(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I666)
   : ta0_(proj), ta1_(I666) { }
};

class Task510 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task510(std::shared_ptr<TATensor<double,2>> I666, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I667)
   : ta0_(I666), ta1_(Gamma38), ta2_(I667) { }
};

class Task511 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task511(std::shared_ptr<TATensor<double,4>> I667, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I667), ta1_(t2) { }
};

class Task512 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task512(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I670)
   : ta0_(proj), ta1_(I670) { }
};

class Task513 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task513(std::shared_ptr<TATensor<double,2>> I670, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I670), ta1_(Gamma60), ta2_(t2) { }
};

class Task514 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task514(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task515 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task515(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I672)
   : ta0_(proj), ta1_(I672) { }
};

class Task516 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task516(std::shared_ptr<TATensor<double,4>> I672, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I672), ta1_(Gamma92), ta2_(t2) { }
};

class Task517 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task517(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I674)
   : ta0_(proj), ta1_(I674) { }
};

class Task518 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task518(std::shared_ptr<TATensor<double,4>> I674, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I674), ta1_(Gamma6), ta2_(t2) { }
};

class Task519 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task519(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I676)
   : ta0_(proj), ta1_(I676) { }
};

class Task520 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task520(std::shared_ptr<TATensor<double,4>> I676, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> I677)
   : ta0_(I676), ta1_(Gamma16), ta2_(I677) { }
};

class Task521 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task521(std::shared_ptr<TATensor<double,4>> I677, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I677), ta1_(t2) { }
};

class Task522 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task522(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I680)
   : ta0_(proj), ta1_(I680) { }
};

class Task523 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task523(std::shared_ptr<TATensor<double,4>> I680, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I680), ta1_(Gamma32), ta2_(t2) { }
};

class Task524 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task524(std::shared_ptr<TATensor<double,4>> I680, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I680), ta1_(Gamma35), ta2_(t2) { }
};

class Task525 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task525(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I684)
   : ta0_(proj), ta1_(I684) { }
};

class Task526 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task526(std::shared_ptr<TATensor<double,4>> I684, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I685)
   : ta0_(I684), ta1_(Gamma35), ta2_(I685) { }
};

class Task527 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task527(std::shared_ptr<TATensor<double,4>> I685, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I685), ta1_(t2) { }
};

class Task528 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task528(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I688)
   : ta0_(proj), ta1_(I688) { }
};

class Task529 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task529(std::shared_ptr<TATensor<double,4>> I688, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I688), ta1_(Gamma59), ta2_(t2) { }
};

class Task530 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task530(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I690)
   : ta0_(proj), ta1_(I690) { }
};

class Task531 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task531(std::shared_ptr<TATensor<double,4>> I690, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I690), ta1_(t2) { }
};

class Task532 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task532(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I692)
   : ta0_(proj), ta1_(I692) { }
};

class Task533 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task533(std::shared_ptr<TATensor<double,4>> I692, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I693)
   : ta0_(I692), ta1_(Gamma38), ta2_(I693) { }
};

class Task534 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task534(std::shared_ptr<TATensor<double,4>> I693, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I693), ta1_(t2) { }
};

class Task535 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task535(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I696)
   : ta0_(proj), ta1_(I696) { }
};

class Task536 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task536(std::shared_ptr<TATensor<double,4>> I696, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I696), ta1_(Gamma60), ta2_(t2) { }
};

class Task537 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task537(std::shared_ptr<TATensor<double,1>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task538 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,1>> ta1_;
    void compute_() override;
  public:
    Task538(std::shared_ptr<TATensor<double,1>> proj, std::shared_ptr<TATensor<double,1>> I698)
   : ta0_(proj), ta1_(I698) { }
};

class Task539 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task539(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma248, std::shared_ptr<TATensor<double,4>> I699)
   : ta0_(I698), ta1_(Gamma248), ta2_(I699) { }
};

class Task540 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task540(std::shared_ptr<TATensor<double,4>> I699, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I699), ta1_(t2) { }
};

class Task541 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task541(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma249, std::shared_ptr<TATensor<double,4>> I702)
   : ta0_(I698), ta1_(Gamma249), ta2_(I702) { }
};

class Task542 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task542(std::shared_ptr<TATensor<double,4>> I702, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I703)
   : ta0_(I702), ta1_(t2), ta2_(I703) { }
};

class Task543 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task543(std::shared_ptr<TATensor<double,4>> I703, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I703), ta1_(f1), ta2_(t2) { }
};

class Task544 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task544(std::shared_ptr<TATensor<double,4>> I702, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I702), ta1_(t2), e0_(e) { }
};

class Task545 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task545(std::shared_ptr<TATensor<double,4>> I702, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I702), ta1_(v2), ta2_(t2) { }
};

class Task546 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task546(std::shared_ptr<TATensor<double,4>> I702, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I702), ta1_(v2), ta2_(t2) { }
};

class Task547 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task547(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma250, std::shared_ptr<TATensor<double,6>> I706)
   : ta0_(I698), ta1_(Gamma250), ta2_(I706) { }
};

class Task548 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task548(std::shared_ptr<TATensor<double,6>> I706, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I707)
   : ta0_(I706), ta1_(t2), ta2_(I707) { }
};

class Task549 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task549(std::shared_ptr<TATensor<double,4>> I707, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I707), ta1_(f1), ta2_(t2) { }
};


}
}
}
#endif
#endif

