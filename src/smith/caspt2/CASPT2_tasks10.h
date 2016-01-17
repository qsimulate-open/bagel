//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks10.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS10_H
#define __SRC_SMITH_CASPT2_TASKS10_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task450 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task450(std::shared_ptr<TATensor<double,2>> I566, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I567)
   : ta0_(I566), ta1_(Gamma16), ta2_(I567) { }
};

class Task451 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task451(std::shared_ptr<TATensor<double,2>> I567, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I567), ta1_(t2) { }
};

class Task452 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task452(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I572)
   : ta0_(proj), ta1_(I572) { }
};

class Task453 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task453(std::shared_ptr<TATensor<double,2>> I572, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I573)
   : ta0_(I572), ta1_(t2), ta2_(I573) { }
};

class Task454 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task454(std::shared_ptr<TATensor<double,2>> I573, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I574)
   : ta0_(I573), ta1_(Gamma38), ta2_(I574) { }
};

class Task455 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task455(std::shared_ptr<TATensor<double,4>> I574, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I574), ta1_(t2) { }
};

class Task456 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task456(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I581)
   : ta0_(proj), ta1_(I581) { }
};

class Task457 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,0>> ta2_;
    void compute_() override;
  public:
    Task457(std::shared_ptr<TATensor<double,2>> I581, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,0>> I582)
   : ta0_(I581), ta1_(Gamma38), ta2_(I582) { }
};

class Task458 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task458(std::shared_ptr<TATensor<double,0>> I582, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I582), ta1_(t2) { }
};

class Task459 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task459(std::shared_ptr<TATensor<double,0>> I582, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I582), ta1_(t2) { }
};

class Task460 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task460(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I587)
   : ta0_(proj), ta1_(I587) { }
};

class Task461 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task461(std::shared_ptr<TATensor<double,2>> I587, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I587), ta1_(t2) { }
};

class Task462 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task462(std::shared_ptr<TATensor<double,2>> I587, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I587), ta1_(t2) { }
};

class Task463 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task463(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I591)
   : ta0_(proj), ta1_(I591) { }
};

class Task464 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task464(std::shared_ptr<TATensor<double,2>> I591, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I591), ta1_(t2) { }
};

class Task465 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task465(std::shared_ptr<TATensor<double,2>> I591, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I591), ta1_(t2) { }
};

class Task466 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task466(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I595)
   : ta0_(proj), ta1_(I595) { }
};

class Task467 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task467(std::shared_ptr<TATensor<double,2>> I595, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I596)
   : ta0_(I595), ta1_(Gamma38), ta2_(I596) { }
};

class Task468 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task468(std::shared_ptr<TATensor<double,2>> I596, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I596), ta1_(t2) { }
};

class Task469 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task469(std::shared_ptr<TATensor<double,2>> I596, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I596), ta1_(t2) { }
};

class Task470 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task470(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I601)
   : ta0_(proj), ta1_(I601) { }
};

class Task471 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task471(std::shared_ptr<TATensor<double,2>> I601, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I602)
   : ta0_(I601), ta1_(t2), ta2_(I602) { }
};

class Task472 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task472(std::shared_ptr<TATensor<double,4>> I602, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I603)
   : ta0_(I602), ta1_(Gamma35), ta2_(I603) { }
};

class Task473 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task473(std::shared_ptr<TATensor<double,4>> I603, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I603), ta1_(t2) { }
};

class Task474 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task474(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I604)
   : ta0_(proj), ta1_(I604) { }
};

class Task475 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task475(std::shared_ptr<TATensor<double,2>> I604, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I605)
   : ta0_(I604), ta1_(t2), ta2_(I605) { }
};

class Task476 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task476(std::shared_ptr<TATensor<double,4>> I605, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I605), ta1_(Gamma32), ta2_(t2) { }
};

class Task477 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task477(std::shared_ptr<TATensor<double,4>> I605, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I605), ta1_(Gamma35), ta2_(t2) { }
};

class Task478 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task478(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I613)
   : ta0_(proj), ta1_(I613) { }
};

class Task479 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task479(std::shared_ptr<TATensor<double,2>> I613, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I614)
   : ta0_(I613), ta1_(t2), ta2_(I614) { }
};

class Task480 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task480(std::shared_ptr<TATensor<double,2>> I614, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I614), ta1_(Gamma60), ta2_(t2) { }
};

class Task481 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task481(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I616)
   : ta0_(proj), ta1_(I616) { }
};

class Task482 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task482(std::shared_ptr<TATensor<double,2>> I616, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I617)
   : ta0_(I616), ta1_(t2), ta2_(I617) { }
};

class Task483 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task483(std::shared_ptr<TATensor<double,2>> I617, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I617), ta1_(Gamma60), ta2_(t2) { }
};

class Task484 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task484(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I619)
   : ta0_(proj), ta1_(I619) { }
};

class Task485 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task485(std::shared_ptr<TATensor<double,2>> I619, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I620)
   : ta0_(I619), ta1_(Gamma38), ta2_(I620) { }
};

class Task486 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task486(std::shared_ptr<TATensor<double,2>> I620, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I620), ta1_(t2) { }
};

class Task487 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task487(std::shared_ptr<TATensor<double,2>> I620, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I620), ta1_(t2) { }
};

class Task488 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task488(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I631)
   : ta0_(proj), ta1_(I631) { }
};

class Task489 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task489(std::shared_ptr<TATensor<double,2>> I631, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I632)
   : ta0_(I631), ta1_(t2), ta2_(I632) { }
};

class Task490 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task490(std::shared_ptr<TATensor<double,4>> I632, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I632), ta1_(Gamma38), ta2_(t2) { }
};

class Task491 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task491(std::shared_ptr<TATensor<double,2>> I631, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I635)
   : ta0_(I631), ta1_(t2), ta2_(I635) { }
};

class Task492 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task492(std::shared_ptr<TATensor<double,4>> I635, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I635), ta1_(Gamma38), ta2_(t2) { }
};

class Task493 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task493(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I637)
   : ta0_(proj), ta1_(I637) { }
};

class Task494 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task494(std::shared_ptr<TATensor<double,2>> I637, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I638)
   : ta0_(I637), ta1_(t2), ta2_(I638) { }
};

class Task495 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task495(std::shared_ptr<TATensor<double,4>> I638, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I638), ta1_(Gamma38), ta2_(t2) { }
};

class Task496 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task496(std::shared_ptr<TATensor<double,2>> I637, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I641)
   : ta0_(I637), ta1_(t2), ta2_(I641) { }
};

class Task497 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task497(std::shared_ptr<TATensor<double,4>> I641, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I641), ta1_(Gamma38), ta2_(t2) { }
};

class Task498 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task498(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I643)
   : ta0_(proj), ta1_(I643) { }
};

class Task499 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task499(std::shared_ptr<TATensor<double,2>> I643, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I644)
   : ta0_(I643), ta1_(t2), ta2_(I644) { }
};


}
}
}
#endif
#endif

