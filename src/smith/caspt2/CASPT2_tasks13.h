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
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task600(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,3>> Gamma264, std::shared_ptr<TATensor<double,2>> I759)
   : ta0_(I698), ta1_(Gamma264), ta2_(I759) { }
};

class Task601 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task601(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I760)
   : ta0_(I759), ta1_(t2), ta2_(I760) { }
};

class Task602 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task602(std::shared_ptr<TATensor<double,4>> I760, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I760), ta1_(f1), ta2_(t2) { }
};

class Task603 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task603(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I764)
   : ta0_(I759), ta1_(t2), ta2_(I764) { }
};

class Task604 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task604(std::shared_ptr<TATensor<double,4>> I764, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I764), ta1_(f1), ta2_(t2) { }
};

class Task605 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task605(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I768)
   : ta0_(I759), ta1_(t2), ta2_(I768) { }
};

class Task606 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task606(std::shared_ptr<TATensor<double,4>> I768, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I768), ta1_(f1), ta2_(t2) { }
};

class Task607 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task607(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I772)
   : ta0_(I759), ta1_(t2), ta2_(I772) { }
};

class Task608 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task608(std::shared_ptr<TATensor<double,4>> I772, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I772), ta1_(f1), ta2_(t2) { }
};

class Task609 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task609(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I776)
   : ta0_(I759), ta1_(t2), ta2_(I776) { }
};

class Task610 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task610(std::shared_ptr<TATensor<double,4>> I776, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I776), ta1_(f1), ta2_(t2) { }
};

class Task611 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task611(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I780)
   : ta0_(I759), ta1_(t2), ta2_(I780) { }
};

class Task612 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task612(std::shared_ptr<TATensor<double,4>> I780, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I780), ta1_(f1), ta2_(t2) { }
};

class Task613 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task613(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I800)
   : ta0_(I759), ta1_(f1), ta2_(I800) { }
};

class Task614 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task614(std::shared_ptr<TATensor<double,2>> I800, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I800), ta1_(t2) { }
};

class Task615 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task615(std::shared_ptr<TATensor<double,2>> I800, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I800), ta1_(t2) { }
};

class Task616 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task616(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I943)
   : ta0_(I759), ta1_(f1), ta2_(I943) { }
};

class Task617 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task617(std::shared_ptr<TATensor<double,2>> I943, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I943), ta1_(t2) { }
};

class Task618 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task618(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I947)
   : ta0_(I759), ta1_(f1), ta2_(I947) { }
};

class Task619 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task619(std::shared_ptr<TATensor<double,2>> I947, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I947), ta1_(t2) { }
};

class Task620 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task620(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I759), ta1_(t2), e0_(e) { }
};

class Task621 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task621(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I759), ta1_(t2), e0_(e) { }
};

class Task622 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task622(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I759), ta1_(v2), ta2_(t2) { }
};

class Task623 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task623(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I759), ta1_(v2), ta2_(t2) { }
};

class Task624 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task624(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I759), ta1_(v2), ta2_(t2) { }
};

class Task625 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task625(std::shared_ptr<TATensor<double,2>> I759, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I759), ta1_(v2), ta2_(t2) { }
};

class Task626 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task626(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma270, std::shared_ptr<TATensor<double,4>> I783)
   : ta0_(I698), ta1_(Gamma270), ta2_(I783) { }
};

class Task627 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task627(std::shared_ptr<TATensor<double,4>> I783, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I784)
   : ta0_(I783), ta1_(t2), ta2_(I784) { }
};

class Task628 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task628(std::shared_ptr<TATensor<double,4>> I784, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I784), ta1_(f1), ta2_(t2) { }
};

class Task629 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task629(std::shared_ptr<TATensor<double,4>> I783, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I783), ta1_(v2), ta2_(t2) { }
};

class Task630 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task630(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma276, std::shared_ptr<TATensor<double,6>> I807)
   : ta0_(I698), ta1_(Gamma276), ta2_(I807) { }
};

class Task631 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task631(std::shared_ptr<TATensor<double,6>> I807, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I808)
   : ta0_(I807), ta1_(t2), ta2_(I808) { }
};

class Task632 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task632(std::shared_ptr<TATensor<double,4>> I808, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I808), ta1_(f1), ta2_(t2) { }
};

class Task633 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task633(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma277, std::shared_ptr<TATensor<double,4>> I811)
   : ta0_(I698), ta1_(Gamma277), ta2_(I811) { }
};

class Task634 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task634(std::shared_ptr<TATensor<double,4>> I811, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I812)
   : ta0_(I811), ta1_(t2), ta2_(I812) { }
};

class Task635 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task635(std::shared_ptr<TATensor<double,4>> I812, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I812), ta1_(t2), ta2_(f1) { }
};

class Task636 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task636(std::shared_ptr<TATensor<double,4>> I811, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I811), ta1_(v2), ta2_(t2) { }
};

class Task637 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task637(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma279, std::shared_ptr<TATensor<double,4>> I819)
   : ta0_(I698), ta1_(Gamma279), ta2_(I819) { }
};

class Task638 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task638(std::shared_ptr<TATensor<double,4>> I819, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I819), ta1_(t2) { }
};

class Task639 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task639(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma280, std::shared_ptr<TATensor<double,4>> I822)
   : ta0_(I698), ta1_(Gamma280), ta2_(I822) { }
};

class Task640 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task640(std::shared_ptr<TATensor<double,4>> I822, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I823)
   : ta0_(I822), ta1_(t2), ta2_(I823) { }
};

class Task641 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task641(std::shared_ptr<TATensor<double,4>> I823, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I823), ta1_(f1), ta2_(t2) { }
};

class Task642 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task642(std::shared_ptr<TATensor<double,4>> I822, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I827)
   : ta0_(I822), ta1_(t2), ta2_(I827) { }
};

class Task643 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task643(std::shared_ptr<TATensor<double,4>> I827, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I827), ta1_(f1), ta2_(t2) { }
};

class Task644 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task644(std::shared_ptr<TATensor<double,4>> I822, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I858)
   : ta0_(I822), ta1_(t2), ta2_(I858) { }
};

class Task645 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task645(std::shared_ptr<TATensor<double,4>> I858, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I858), ta1_(t2), ta2_(f1) { }
};

class Task646 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task646(std::shared_ptr<TATensor<double,4>> I822, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I985)
   : ta0_(I822), ta1_(t2), ta2_(I985) { }
};

class Task647 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task647(std::shared_ptr<TATensor<double,4>> I985, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I985), ta1_(t2), e0_(e) { }
};

class Task648 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task648(std::shared_ptr<TATensor<double,4>> I985, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I985), ta1_(f1), ta2_(t2) { }
};

class Task649 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task649(std::shared_ptr<TATensor<double,4>> I822, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I822), ta1_(v2), ta2_(t2) { }
};


}
}
}
#endif
#endif

