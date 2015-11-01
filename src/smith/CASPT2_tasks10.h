//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks10.h
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
    void compute_() override;
  public:
    Task450(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I585)
   : ta0_(proj), ta1_(I585) { }
};

class Task451 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task451(std::shared_ptr<TATensor<double,2>> I585, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I586)
   : ta0_(I585), ta1_(t2), ta2_(I586) { }
};

class Task452 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task452(std::shared_ptr<TATensor<double,4>> I586, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I586), ta1_(Gamma35), ta2_(t2) { }
};

class Task453 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task453(std::shared_ptr<TATensor<double,2>> I585, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I595)
   : ta0_(I585), ta1_(t2), ta2_(I595) { }
};

class Task454 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task454(std::shared_ptr<TATensor<double,4>> I595, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I595), ta1_(Gamma35), ta2_(t2) { }
};

class Task455 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task455(std::shared_ptr<TATensor<double,2>> I585, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I732)
   : ta0_(I585), ta1_(t2), ta2_(I732) { }
};

class Task456 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task456(std::shared_ptr<TATensor<double,4>> I732, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I732), ta1_(Gamma60), ta2_(t2) { }
};

class Task457 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task457(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I597)
   : ta0_(proj), ta1_(I597) { }
};

class Task458 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task458(std::shared_ptr<TATensor<double,2>> I597, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I598)
   : ta0_(I597), ta1_(t2), ta2_(I598) { }
};

class Task459 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task459(std::shared_ptr<TATensor<double,4>> I598, std::shared_ptr<TATensor<double,6>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I598), ta1_(Gamma51), ta2_(t2) { }
};

class Task460 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task460(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I621)
   : ta0_(proj), ta1_(I621) { }
};

class Task461 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task461(std::shared_ptr<TATensor<double,2>> I621, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I622)
   : ta0_(I621), ta1_(t2), ta2_(I622) { }
};

class Task462 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task462(std::shared_ptr<TATensor<double,4>> I622, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I622), ta1_(Gamma59), ta2_(t2) { }
};

class Task463 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task463(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I633)
   : ta0_(proj), ta1_(I633) { }
};

class Task464 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task464(std::shared_ptr<TATensor<double,2>> I633, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I634)
   : ta0_(I633), ta1_(Gamma16), ta2_(I634) { }
};

class Task465 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task465(std::shared_ptr<TATensor<double,2>> I634, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I634), ta1_(t2) { }
};

class Task466 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task466(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I636)
   : ta0_(proj), ta1_(I636) { }
};

class Task467 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task467(std::shared_ptr<TATensor<double,2>> I636, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I637)
   : ta0_(I636), ta1_(Gamma16), ta2_(I637) { }
};

class Task468 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task468(std::shared_ptr<TATensor<double,2>> I637, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I637), ta1_(t2) { }
};

class Task469 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task469(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I642)
   : ta0_(proj), ta1_(I642) { }
};

class Task470 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task470(std::shared_ptr<TATensor<double,2>> I642, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I643)
   : ta0_(I642), ta1_(t2), ta2_(I643) { }
};

class Task471 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task471(std::shared_ptr<TATensor<double,2>> I643, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I644)
   : ta0_(I643), ta1_(Gamma38), ta2_(I644) { }
};

class Task472 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task472(std::shared_ptr<TATensor<double,4>> I644, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I644), ta1_(t2) { }
};

class Task473 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task473(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I651)
   : ta0_(proj), ta1_(I651) { }
};

class Task474 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,0>> ta2_;
    void compute_() override;
  public:
    Task474(std::shared_ptr<TATensor<double,2>> I651, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,0>> I652)
   : ta0_(I651), ta1_(Gamma38), ta2_(I652) { }
};

class Task475 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task475(std::shared_ptr<TATensor<double,0>> I652, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I652), ta1_(t2) { }
};

class Task476 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task476(std::shared_ptr<TATensor<double,0>> I652, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I652), ta1_(t2) { }
};

class Task477 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task477(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I657)
   : ta0_(proj), ta1_(I657) { }
};

class Task478 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task478(std::shared_ptr<TATensor<double,2>> I657, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I657), ta1_(t2) { }
};

class Task479 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task479(std::shared_ptr<TATensor<double,2>> I657, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I657), ta1_(t2) { }
};

class Task480 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task480(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I661)
   : ta0_(proj), ta1_(I661) { }
};

class Task481 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task481(std::shared_ptr<TATensor<double,2>> I661, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I661), ta1_(t2) { }
};

class Task482 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task482(std::shared_ptr<TATensor<double,2>> I661, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I661), ta1_(t2) { }
};

class Task483 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task483(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I665)
   : ta0_(proj), ta1_(I665) { }
};

class Task484 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task484(std::shared_ptr<TATensor<double,2>> I665, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I666)
   : ta0_(I665), ta1_(Gamma38), ta2_(I666) { }
};

class Task485 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task485(std::shared_ptr<TATensor<double,2>> I666, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I666), ta1_(t2) { }
};

class Task486 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task486(std::shared_ptr<TATensor<double,2>> I666, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I666), ta1_(t2) { }
};

class Task487 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task487(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I671)
   : ta0_(proj), ta1_(I671) { }
};

class Task488 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task488(std::shared_ptr<TATensor<double,2>> I671, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I672)
   : ta0_(I671), ta1_(t2), ta2_(I672) { }
};

class Task489 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task489(std::shared_ptr<TATensor<double,4>> I672, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I673)
   : ta0_(I672), ta1_(Gamma35), ta2_(I673) { }
};

class Task490 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task490(std::shared_ptr<TATensor<double,4>> I673, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I673), ta1_(t2) { }
};

class Task491 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task491(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I674)
   : ta0_(proj), ta1_(I674) { }
};

class Task492 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task492(std::shared_ptr<TATensor<double,2>> I674, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I675)
   : ta0_(I674), ta1_(t2), ta2_(I675) { }
};

class Task493 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task493(std::shared_ptr<TATensor<double,4>> I675, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I675), ta1_(Gamma32), ta2_(t2) { }
};

class Task494 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task494(std::shared_ptr<TATensor<double,4>> I675, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I675), ta1_(Gamma35), ta2_(t2) { }
};

class Task495 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task495(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I683)
   : ta0_(proj), ta1_(I683) { }
};

class Task496 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task496(std::shared_ptr<TATensor<double,2>> I683, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I684)
   : ta0_(I683), ta1_(t2), ta2_(I684) { }
};

class Task497 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task497(std::shared_ptr<TATensor<double,2>> I684, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I684), ta1_(Gamma60), ta2_(t2) { }
};

class Task498 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task498(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I686)
   : ta0_(proj), ta1_(I686) { }
};

class Task499 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task499(std::shared_ptr<TATensor<double,2>> I686, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I687)
   : ta0_(I686), ta1_(t2), ta2_(I687) { }
};


}
}
}
#endif
#endif

