//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks6.h
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

class Task250 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task250(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(v2), ta2_(t2) { }
};

class Task251 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task251(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(v2), ta2_(t2) { }
};

class Task252 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task252(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(v2), ta2_(t2) { }
};

class Task253 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task253(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(v2), ta2_(t2) { }
};

class Task254 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task254(std::shared_ptr<TATensor<double,4>> I350, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I350), ta1_(v2), ta2_(t2) { }
};

class Task255 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task255(std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I353)
   : ta1_(Gamma29), ta2_(I353) { }
};

class Task256 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task256(std::shared_ptr<TATensor<double,4>> I353, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I353), ta1_(v2), ta2_(t2) { }
};

class Task257 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task257(std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> I356)
   : ta1_(Gamma32), ta2_(I356) { }
};

class Task258 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task258(std::shared_ptr<TATensor<double,4>> I356, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I356), ta1_(v2), ta2_(t2) { }
};

class Task259 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task259(std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> I365)
   : ta1_(Gamma7), ta2_(I365) { }
};

class Task260 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task260(std::shared_ptr<TATensor<double,4>> I365, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I365), ta1_(v2), ta2_(t2) { }
};

class Task261 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task261(std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,6>> I374)
   : ta1_(Gamma59), ta2_(I374) { }
};

class Task262 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task262(std::shared_ptr<TATensor<double,6>> I374, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I374), ta1_(v2), ta2_(t2) { }
};

class Task263 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task263(std::shared_ptr<TATensor<double,6>> Gamma57, std::shared_ptr<TATensor<double,6>> I377)
   : ta1_(Gamma57), ta2_(I377) { }
};

class Task264 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task264(std::shared_ptr<TATensor<double,6>> I377, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I377), ta1_(v2), ta2_(t2) { }
};

class Task265 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task265(std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta1_(v2), ta2_(t2) { }
};

class Task266 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task266(std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta1_(v2), ta2_(t2) { }
};

class Task267 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task267(std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I384)
   : ta1_(Gamma38), ta2_(I384) { }
};

class Task268 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task268(std::shared_ptr<TATensor<double,2>> I384, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I384), ta1_(v2), ta2_(t2) { }
};

class Task269 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task269(std::shared_ptr<TATensor<double,2>> I384, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I384), ta1_(v2), ta2_(t2) { }
};

class Task270 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task270(std::shared_ptr<TATensor<double,2>> I384, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I384), ta1_(h1), ta2_(t2) { }
};

class Task271 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task271(std::shared_ptr<TATensor<double,2>> I384, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I384), ta1_(h1), ta2_(t2) { }
};

class Task272 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task272(std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> I390)
   : ta1_(Gamma60), ta2_(I390) { }
};

class Task273 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task273(std::shared_ptr<TATensor<double,4>> I390, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I390), ta1_(v2), ta2_(t2) { }
};

class Task274 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task274(std::shared_ptr<TATensor<double,4>> I390, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I390), ta1_(h1), ta2_(t2) { }
};

class Task275 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task275(std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I393)
   : ta1_(h1), ta2_(I393) { }
};

class Task276 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task276(std::shared_ptr<TATensor<double,2>> I393, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I393), ta1_(Gamma7), ta2_(t2) { }
};

class Task277 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task277(std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> I405)
   : ta1_(Gamma92), ta2_(I405) { }
};

class Task278 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task278(std::shared_ptr<TATensor<double,4>> I405, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I405), ta1_(t2) { }
};

class Task279 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task279(std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I408)
   : ta1_(t2), ta2_(I408) { }
};

class Task280 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task280(std::shared_ptr<TATensor<double,4>> I408, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I408), ta1_(Gamma6), ta2_(t2) { }
};

class Task281 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task281(std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I411)
   : ta1_(Gamma16), ta2_(I411) { }
};

class Task282 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task282(std::shared_ptr<TATensor<double,2>> I411, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I411), ta1_(t2) { }
};

class Task283 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task283(std::shared_ptr<TATensor<double,2>> I411, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I411), ta1_(t2) { }
};

class Task284 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task284(std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> I417)
   : ta1_(Gamma32), ta2_(I417) { }
};

class Task285 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task285(std::shared_ptr<TATensor<double,4>> I417, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I417), ta1_(t2) { }
};

class Task286 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task286(std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I420)
   : ta1_(Gamma35), ta2_(I420) { }
};

class Task287 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task287(std::shared_ptr<TATensor<double,4>> I420, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I420), ta1_(t2) { }
};

class Task288 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task288(std::shared_ptr<TATensor<double,4>> I420, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I420), ta1_(t2) { }
};

class Task289 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task289(std::shared_ptr<TATensor<double,4>> I420, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I420), ta1_(t2) { }
};

class Task290 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task290(std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,6>> I429)
   : ta1_(Gamma59), ta2_(I429) { }
};

class Task291 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task291(std::shared_ptr<TATensor<double,6>> I429, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I429), ta1_(t2) { }
};

class Task292 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task292(std::shared_ptr<TATensor<double,4>> t2)
   : ta1_(t2) { }
};

class Task293 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task293(std::shared_ptr<TATensor<double,4>> t2)
   : ta1_(t2) { }
};

class Task294 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task294(std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I436)
   : ta1_(Gamma38), ta2_(I436) { }
};

class Task295 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task295(std::shared_ptr<TATensor<double,2>> I436, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I436), ta1_(t2) { }
};

class Task296 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task296(std::shared_ptr<TATensor<double,2>> I436, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I436), ta1_(t2) { }
};

class Task297 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task297(std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> I442)
   : ta1_(Gamma60), ta2_(I442) { }
};

class Task298 : public AccTask {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task298(std::shared_ptr<TATensor<double,4>> I442, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I442), ta1_(t2) { }
};

class Task299 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task299(std::shared_ptr<TATensor<double,2>> t, const bool reset) : tensor_(t), reset_(reset) { }
};


}
}
}
#endif
#endif

