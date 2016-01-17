//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks6.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS6_H
#define __SRC_SMITH_CASPT2_TASKS6_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task250 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task250(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I314)
   : ta0_(proj), ta1_(I314) { }
};

class Task251 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task251(std::shared_ptr<TATensor<double,4>> I314, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I314), ta1_(Gamma59), ta2_(v2) { }
};

class Task252 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task252(std::shared_ptr<TATensor<double,4>> I314, std::shared_ptr<TATensor<double,6>> Gamma57, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I314), ta1_(Gamma57), ta2_(v2) { }
};

class Task253 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task253(std::shared_ptr<TATensor<double,4>> I314, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I314), ta1_(Gamma60), ta2_(h1) { }
};

class Task254 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task254(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I318)
   : ta0_(proj), ta1_(I318) { }
};

class Task255 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task255(std::shared_ptr<TATensor<double,4>> I318, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I318), ta1_(v2) { }
};

class Task256 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task256(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I320)
   : ta0_(proj), ta1_(I320) { }
};

class Task257 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task257(std::shared_ptr<TATensor<double,4>> I320, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I321)
   : ta0_(I320), ta1_(Gamma38), ta2_(I321) { }
};

class Task258 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task258(std::shared_ptr<TATensor<double,4>> I321, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I321), ta1_(v2) { }
};

class Task259 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task259(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I324)
   : ta0_(proj), ta1_(I324) { }
};

class Task260 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task260(std::shared_ptr<TATensor<double,4>> I324, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I324), ta1_(Gamma60), ta2_(v2) { }
};

class Task261 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task261(std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> I335)
   : ta1_(Gamma92), ta2_(I335) { }
};

class Task262 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task262(std::shared_ptr<TATensor<double,4>> I335, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I335), ta1_(t2) { }
};

class Task263 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task263(std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I338)
   : ta1_(t2), ta2_(I338) { }
};

class Task264 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task264(std::shared_ptr<TATensor<double,4>> I338, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I338), ta1_(Gamma6), ta2_(t2) { }
};

class Task265 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task265(std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I341)
   : ta1_(Gamma16), ta2_(I341) { }
};

class Task266 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task266(std::shared_ptr<TATensor<double,2>> I341, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I341), ta1_(t2) { }
};

class Task267 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task267(std::shared_ptr<TATensor<double,2>> I341, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I341), ta1_(t2) { }
};

class Task268 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task268(std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> I347)
   : ta1_(Gamma32), ta2_(I347) { }
};

class Task269 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task269(std::shared_ptr<TATensor<double,4>> I347, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I347), ta1_(t2) { }
};

class Task270 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task270(std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I350)
   : ta1_(Gamma35), ta2_(I350) { }
};

class Task271 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task271(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(t2) { }
};

class Task272 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task272(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(t2) { }
};

class Task273 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task273(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(t2) { }
};

class Task274 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task274(std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,6>> I359)
   : ta1_(Gamma59), ta2_(I359) { }
};

class Task275 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task275(std::shared_ptr<TATensor<double,6>> I359, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I359), ta1_(t2) { }
};

class Task276 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task276(std::shared_ptr<TATensor<double,4>> t2)
   : ta1_(t2) { }
};

class Task277 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task277(std::shared_ptr<TATensor<double,4>> t2)
   : ta1_(t2) { }
};

class Task278 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task278(std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I366)
   : ta1_(Gamma38), ta2_(I366) { }
};

class Task279 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task279(std::shared_ptr<TATensor<double,2>> I366, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I366), ta1_(t2) { }
};

class Task280 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task280(std::shared_ptr<TATensor<double,2>> I366, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I366), ta1_(t2) { }
};

class Task281 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task281(std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> I372)
   : ta1_(Gamma60), ta2_(I372) { }
};

class Task282 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task282(std::shared_ptr<TATensor<double,4>> I372, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I372), ta1_(t2) { }
};

class Task283 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task283(std::shared_ptr<TATensor<double,2>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task284 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task284(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I374)
   : ta0_(proj), ta1_(I374) { }
};

class Task285 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task285(std::shared_ptr<TATensor<double,2>> I374, std::shared_ptr<TATensor<double,6>> Gamma138, std::shared_ptr<TATensor<double,4>> I375)
   : ta0_(I374), ta1_(Gamma138), ta2_(I375) { }
};

class Task286 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task286(std::shared_ptr<TATensor<double,4>> I375, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I375), ta1_(t2) { }
};

class Task287 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task287(std::shared_ptr<TATensor<double,2>> I374, std::shared_ptr<TATensor<double,6>> Gamma169, std::shared_ptr<TATensor<double,4>> I468)
   : ta0_(I374), ta1_(Gamma169), ta2_(I468) { }
};

class Task288 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task288(std::shared_ptr<TATensor<double,4>> I468, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I468), ta1_(t2) { }
};

class Task289 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task289(std::shared_ptr<TATensor<double,2>> I374, std::shared_ptr<TATensor<double,6>> Gamma172, std::shared_ptr<TATensor<double,4>> I477)
   : ta0_(I374), ta1_(Gamma172), ta2_(I477) { }
};

class Task290 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task290(std::shared_ptr<TATensor<double,4>> I477, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I477), ta1_(t2) { }
};

class Task291 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task291(std::shared_ptr<TATensor<double,4>> I477, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I477), ta1_(t2) { }
};

class Task292 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task292(std::shared_ptr<TATensor<double,4>> I477, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I477), ta1_(t2) { }
};

class Task293 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task293(std::shared_ptr<TATensor<double,2>> I374, std::shared_ptr<TATensor<double,6>> Gamma230, std::shared_ptr<TATensor<double,4>> I659)
   : ta0_(I374), ta1_(Gamma230), ta2_(I659) { }
};

class Task294 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task294(std::shared_ptr<TATensor<double,4>> I659, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I659), ta1_(t2) { }
};

class Task295 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task295(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I377)
   : ta0_(proj), ta1_(I377) { }
};

class Task296 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task296(std::shared_ptr<TATensor<double,2>> I377, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I378)
   : ta0_(I377), ta1_(t2), ta2_(I378) { }
};

class Task297 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task297(std::shared_ptr<TATensor<double,4>> I378, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I378), ta1_(Gamma92), ta2_(t2) { }
};

class Task298 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task298(std::shared_ptr<TATensor<double,2>> I377, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I471)
   : ta0_(I377), ta1_(t2), ta2_(I471) { }
};

class Task299 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task299(std::shared_ptr<TATensor<double,4>> I471, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I471), ta1_(Gamma32), ta2_(t2) { }
};


}
}
}
#endif
#endif

