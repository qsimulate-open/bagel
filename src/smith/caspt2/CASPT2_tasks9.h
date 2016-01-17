//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks9.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS9_H
#define __SRC_SMITH_CASPT2_TASKS9_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task400 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task400(std::shared_ptr<TATensor<double,4>> I462, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I462), ta1_(Gamma29), ta2_(t2) { }
};

class Task401 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task401(std::shared_ptr<TATensor<double,2>> I461, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I465)
   : ta0_(I461), ta1_(t2), ta2_(I465) { }
};

class Task402 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task402(std::shared_ptr<TATensor<double,4>> I465, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I465), ta1_(Gamma7), ta2_(t2) { }
};

class Task403 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task403(std::shared_ptr<TATensor<double,2>> I461, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I504)
   : ta0_(I461), ta1_(t2), ta2_(I504) { }
};

class Task404 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task404(std::shared_ptr<TATensor<double,4>> I504, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I504), ta1_(Gamma7), ta2_(t2) { }
};

class Task405 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task405(std::shared_ptr<TATensor<double,2>> I461, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I507)
   : ta0_(I461), ta1_(t2), ta2_(I507) { }
};

class Task406 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task406(std::shared_ptr<TATensor<double,4>> I507, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I507), ta1_(Gamma7), ta2_(t2) { }
};

class Task407 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task407(std::shared_ptr<TATensor<double,2>> I461, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I656)
   : ta0_(I461), ta1_(t2), ta2_(I656) { }
};

class Task408 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task408(std::shared_ptr<TATensor<double,4>> I656, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I656), ta1_(Gamma60), ta2_(t2) { }
};

class Task409 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task409(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I473)
   : ta0_(proj), ta1_(I473) { }
};

class Task410 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task410(std::shared_ptr<TATensor<double,2>> I473, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I474)
   : ta0_(I473), ta1_(t2), ta2_(I474) { }
};

class Task411 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task411(std::shared_ptr<TATensor<double,4>> I474, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I474), ta1_(Gamma32), ta2_(t2) { }
};

class Task412 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task412(std::shared_ptr<TATensor<double,2>> I473, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I483)
   : ta0_(I473), ta1_(t2), ta2_(I483) { }
};

class Task413 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task413(std::shared_ptr<TATensor<double,4>> I483, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I483), ta1_(Gamma35), ta2_(t2) { }
};

class Task414 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task414(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I488)
   : ta0_(proj), ta1_(I488) { }
};

class Task415 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task415(std::shared_ptr<TATensor<double,2>> I488, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I489)
   : ta0_(I488), ta1_(t2), ta2_(I489) { }
};

class Task416 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task416(std::shared_ptr<TATensor<double,2>> I489, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I489), ta1_(Gamma38), ta2_(t2) { }
};

class Task417 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task417(std::shared_ptr<TATensor<double,2>> I488, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I492)
   : ta0_(I488), ta1_(t2), ta2_(I492) { }
};

class Task418 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task418(std::shared_ptr<TATensor<double,2>> I492, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I492), ta1_(Gamma38), ta2_(t2) { }
};

class Task419 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task419(std::shared_ptr<TATensor<double,2>> I488, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I531)
   : ta0_(I488), ta1_(t2), ta2_(I531) { }
};

class Task420 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task420(std::shared_ptr<TATensor<double,2>> I531, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I531), ta1_(Gamma38), ta2_(t2) { }
};

class Task421 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task421(std::shared_ptr<TATensor<double,2>> I488, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I534)
   : ta0_(I488), ta1_(t2), ta2_(I534) { }
};

class Task422 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task422(std::shared_ptr<TATensor<double,2>> I534, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I534), ta1_(Gamma38), ta2_(t2) { }
};

class Task423 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task423(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I500)
   : ta0_(proj), ta1_(I500) { }
};

class Task424 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task424(std::shared_ptr<TATensor<double,2>> I500, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I501)
   : ta0_(I500), ta1_(t2), ta2_(I501) { }
};

class Task425 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task425(std::shared_ptr<TATensor<double,4>> I501, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I501), ta1_(Gamma6), ta2_(t2) { }
};

class Task426 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task426(std::shared_ptr<TATensor<double,2>> I500, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I653)
   : ta0_(I500), ta1_(t2), ta2_(I653) { }
};

class Task427 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task427(std::shared_ptr<TATensor<double,4>> I653, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I653), ta1_(Gamma59), ta2_(t2) { }
};

class Task428 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task428(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I512)
   : ta0_(proj), ta1_(I512) { }
};

class Task429 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task429(std::shared_ptr<TATensor<double,2>> I512, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I513)
   : ta0_(I512), ta1_(t2), ta2_(I513) { }
};

class Task430 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task430(std::shared_ptr<TATensor<double,4>> I513, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I513), ta1_(Gamma35), ta2_(t2) { }
};

class Task431 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task431(std::shared_ptr<TATensor<double,2>> I512, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I522)
   : ta0_(I512), ta1_(t2), ta2_(I522) { }
};

class Task432 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task432(std::shared_ptr<TATensor<double,4>> I522, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I522), ta1_(Gamma35), ta2_(t2) { }
};

class Task433 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task433(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I515)
   : ta0_(proj), ta1_(I515) { }
};

class Task434 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task434(std::shared_ptr<TATensor<double,2>> I515, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I516)
   : ta0_(I515), ta1_(t2), ta2_(I516) { }
};

class Task435 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task435(std::shared_ptr<TATensor<double,4>> I516, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I516), ta1_(Gamma35), ta2_(t2) { }
};

class Task436 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task436(std::shared_ptr<TATensor<double,2>> I515, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I525)
   : ta0_(I515), ta1_(t2), ta2_(I525) { }
};

class Task437 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task437(std::shared_ptr<TATensor<double,4>> I525, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I525), ta1_(Gamma35), ta2_(t2) { }
};

class Task438 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task438(std::shared_ptr<TATensor<double,2>> I515, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I662)
   : ta0_(I515), ta1_(t2), ta2_(I662) { }
};

class Task439 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task439(std::shared_ptr<TATensor<double,4>> I662, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I662), ta1_(Gamma60), ta2_(t2) { }
};

class Task440 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task440(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I527)
   : ta0_(proj), ta1_(I527) { }
};

class Task441 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task441(std::shared_ptr<TATensor<double,2>> I527, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I528)
   : ta0_(I527), ta1_(t2), ta2_(I528) { }
};

class Task442 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task442(std::shared_ptr<TATensor<double,4>> I528, std::shared_ptr<TATensor<double,6>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I528), ta1_(Gamma51), ta2_(t2) { }
};

class Task443 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task443(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I551)
   : ta0_(proj), ta1_(I551) { }
};

class Task444 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task444(std::shared_ptr<TATensor<double,2>> I551, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I552)
   : ta0_(I551), ta1_(t2), ta2_(I552) { }
};

class Task445 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task445(std::shared_ptr<TATensor<double,4>> I552, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I552), ta1_(Gamma59), ta2_(t2) { }
};

class Task446 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task446(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I563)
   : ta0_(proj), ta1_(I563) { }
};

class Task447 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task447(std::shared_ptr<TATensor<double,2>> I563, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I564)
   : ta0_(I563), ta1_(Gamma16), ta2_(I564) { }
};

class Task448 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task448(std::shared_ptr<TATensor<double,2>> I564, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I564), ta1_(t2) { }
};

class Task449 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task449(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I566)
   : ta0_(proj), ta1_(I566) { }
};


}
}
}
#endif
#endif

