//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks13.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS13_H
#define __SRC_SMITH_CASPT2_TASKS13_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task600 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task600(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma282, std::shared_ptr<TATensor<double,4>> I815)
   : ta0_(I768), ta1_(Gamma282), ta2_(I815) { }
};

class Task601 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task601(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I816)
   : ta0_(I815), ta1_(t2), ta2_(I816) { }
};

class Task602 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task602(std::shared_ptr<TATensor<double,2>> I816, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I816), ta1_(f1), ta2_(t2) { }
};

class Task603 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task603(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I820)
   : ta0_(I815), ta1_(t2), ta2_(I820) { }
};

class Task604 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task604(std::shared_ptr<TATensor<double,2>> I820, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I820), ta1_(f1), ta2_(t2) { }
};

class Task605 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task605(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I858)
   : ta0_(I815), ta1_(t2), ta2_(I858) { }
};

class Task606 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task606(std::shared_ptr<TATensor<double,4>> I858, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I858), ta1_(f1), ta2_(t2) { }
};

class Task607 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task607(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I862)
   : ta0_(I815), ta1_(t2), ta2_(I862) { }
};

class Task608 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task608(std::shared_ptr<TATensor<double,4>> I862, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I862), ta1_(f1), ta2_(t2) { }
};

class Task609 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task609(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I866)
   : ta0_(I815), ta1_(t2), ta2_(I866) { }
};

class Task610 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task610(std::shared_ptr<TATensor<double,4>> I866, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I866), ta1_(f1), ta2_(t2) { }
};

class Task611 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task611(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I815), ta1_(v2), ta2_(t2) { }
};

class Task612 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task612(std::shared_ptr<TATensor<double,4>> I815, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I815), ta1_(h1), ta2_(t2) { }
};

class Task613 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task613(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,3>> Gamma284, std::shared_ptr<TATensor<double,2>> I823)
   : ta0_(I768), ta1_(Gamma284), ta2_(I823) { }
};

class Task614 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task614(std::shared_ptr<TATensor<double,2>> I823, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I823), ta1_(t2) { }
};

class Task615 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task615(std::shared_ptr<TATensor<double,2>> I823, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I823), ta1_(t2) { }
};

class Task616 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task616(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,3>> Gamma286, std::shared_ptr<TATensor<double,2>> I829)
   : ta0_(I768), ta1_(Gamma286), ta2_(I829) { }
};

class Task617 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task617(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I830)
   : ta0_(I829), ta1_(t2), ta2_(I830) { }
};

class Task618 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task618(std::shared_ptr<TATensor<double,4>> I830, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I830), ta1_(f1), ta2_(t2) { }
};

class Task619 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task619(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I834)
   : ta0_(I829), ta1_(t2), ta2_(I834) { }
};

class Task620 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task620(std::shared_ptr<TATensor<double,4>> I834, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I834), ta1_(f1), ta2_(t2) { }
};

class Task621 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task621(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I838)
   : ta0_(I829), ta1_(t2), ta2_(I838) { }
};

class Task622 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task622(std::shared_ptr<TATensor<double,4>> I838, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I838), ta1_(f1), ta2_(t2) { }
};

class Task623 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task623(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I842)
   : ta0_(I829), ta1_(t2), ta2_(I842) { }
};

class Task624 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task624(std::shared_ptr<TATensor<double,4>> I842, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I842), ta1_(f1), ta2_(t2) { }
};

class Task625 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task625(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I846)
   : ta0_(I829), ta1_(t2), ta2_(I846) { }
};

class Task626 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task626(std::shared_ptr<TATensor<double,4>> I846, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I846), ta1_(f1), ta2_(t2) { }
};

class Task627 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task627(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I850)
   : ta0_(I829), ta1_(t2), ta2_(I850) { }
};

class Task628 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task628(std::shared_ptr<TATensor<double,4>> I850, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I850), ta1_(f1), ta2_(t2) { }
};

class Task629 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task629(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I870)
   : ta0_(I829), ta1_(f1), ta2_(I870) { }
};

class Task630 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task630(std::shared_ptr<TATensor<double,2>> I870, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I870), ta1_(t2) { }
};

class Task631 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task631(std::shared_ptr<TATensor<double,2>> I870, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I870), ta1_(t2) { }
};

class Task632 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task632(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I1013)
   : ta0_(I829), ta1_(f1), ta2_(I1013) { }
};

class Task633 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task633(std::shared_ptr<TATensor<double,2>> I1013, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1013), ta1_(t2) { }
};

class Task634 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task634(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I1017)
   : ta0_(I829), ta1_(f1), ta2_(I1017) { }
};

class Task635 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task635(std::shared_ptr<TATensor<double,2>> I1017, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1017), ta1_(t2) { }
};

class Task636 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task636(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I829), ta1_(t2), e0_(e) { }
};

class Task637 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task637(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I829), ta1_(t2), e0_(e) { }
};

class Task638 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task638(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I829), ta1_(v2), ta2_(t2) { }
};

class Task639 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task639(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I829), ta1_(v2), ta2_(t2) { }
};

class Task640 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task640(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I829), ta1_(v2), ta2_(t2) { }
};

class Task641 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task641(std::shared_ptr<TATensor<double,2>> I829, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I829), ta1_(v2), ta2_(t2) { }
};

class Task642 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task642(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma292, std::shared_ptr<TATensor<double,4>> I853)
   : ta0_(I768), ta1_(Gamma292), ta2_(I853) { }
};

class Task643 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task643(std::shared_ptr<TATensor<double,4>> I853, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I854)
   : ta0_(I853), ta1_(t2), ta2_(I854) { }
};

class Task644 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task644(std::shared_ptr<TATensor<double,4>> I854, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I854), ta1_(f1), ta2_(t2) { }
};

class Task645 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task645(std::shared_ptr<TATensor<double,4>> I853, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I853), ta1_(v2), ta2_(t2) { }
};

class Task646 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task646(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma298, std::shared_ptr<TATensor<double,6>> I877)
   : ta0_(I768), ta1_(Gamma298), ta2_(I877) { }
};

class Task647 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task647(std::shared_ptr<TATensor<double,6>> I877, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I878)
   : ta0_(I877), ta1_(t2), ta2_(I878) { }
};

class Task648 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task648(std::shared_ptr<TATensor<double,4>> I878, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I878), ta1_(f1), ta2_(t2) { }
};

class Task649 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task649(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma299, std::shared_ptr<TATensor<double,4>> I881)
   : ta0_(I768), ta1_(Gamma299), ta2_(I881) { }
};


}
}
}
#endif
#endif

