//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks14.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS14_H
#define __SRC_SMITH_CASPT2_TASKS14_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task650 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task650(std::shared_ptr<TATensor<double,4>> I822, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I822), ta1_(v2), ta2_(t2) { }
};

class Task651 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task651(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma282, std::shared_ptr<TATensor<double,4>> I830)
   : ta0_(I698), ta1_(Gamma282), ta2_(I830) { }
};

class Task652 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task652(std::shared_ptr<TATensor<double,4>> I830, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I830), ta1_(t2) { }
};

class Task653 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task653(std::shared_ptr<TATensor<double,4>> I830, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I830), ta1_(t2) { }
};

class Task654 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task654(std::shared_ptr<TATensor<double,4>> I830, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I830), ta1_(t2) { }
};

class Task655 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task655(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma283, std::shared_ptr<TATensor<double,4>> I833)
   : ta0_(I698), ta1_(Gamma283), ta2_(I833) { }
};

class Task656 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task656(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I834)
   : ta0_(I833), ta1_(t2), ta2_(I834) { }
};

class Task657 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task657(std::shared_ptr<TATensor<double,4>> I834, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I834), ta1_(f1), ta2_(t2) { }
};

class Task658 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task658(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I838)
   : ta0_(I833), ta1_(t2), ta2_(I838) { }
};

class Task659 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task659(std::shared_ptr<TATensor<double,4>> I838, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I838), ta1_(f1), ta2_(t2) { }
};

class Task660 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task660(std::shared_ptr<TATensor<double,4>> I838, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I838), ta1_(f1), ta2_(t2) { }
};

class Task661 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task661(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I854)
   : ta0_(I833), ta1_(t2), ta2_(I854) { }
};

class Task662 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task662(std::shared_ptr<TATensor<double,4>> I854, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I854), ta1_(t2), ta2_(f1) { }
};

class Task663 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task663(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I877)
   : ta0_(I833), ta1_(t2), ta2_(I877) { }
};

class Task664 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task664(std::shared_ptr<TATensor<double,4>> I877, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I877), ta1_(f1), ta2_(t2) { }
};

class Task665 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task665(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I881)
   : ta0_(I833), ta1_(t2), ta2_(I881) { }
};

class Task666 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task666(std::shared_ptr<TATensor<double,4>> I881, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I881), ta1_(f1), ta2_(t2) { }
};

class Task667 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task667(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I888)
   : ta0_(I833), ta1_(t2), ta2_(I888) { }
};

class Task668 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task668(std::shared_ptr<TATensor<double,4>> I888, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I888), ta1_(f1), ta2_(t2) { }
};

class Task669 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task669(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I892)
   : ta0_(I833), ta1_(t2), ta2_(I892) { }
};

class Task670 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task670(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I892), ta1_(f1), ta2_(t2) { }
};

class Task671 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task671(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I908)
   : ta0_(I833), ta1_(t2), ta2_(I908) { }
};

class Task672 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task672(std::shared_ptr<TATensor<double,4>> I908, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I908), ta1_(t2), ta2_(f1) { }
};

class Task673 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task673(std::shared_ptr<TATensor<double,4>> I908, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I908), ta1_(t2), ta2_(f1) { }
};

class Task674 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task674(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I981)
   : ta0_(I833), ta1_(t2), ta2_(I981) { }
};

class Task675 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task675(std::shared_ptr<TATensor<double,4>> I981, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I981), ta1_(f1), ta2_(t2) { }
};

class Task676 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task676(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I993)
   : ta0_(I833), ta1_(t2), ta2_(I993) { }
};

class Task677 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task677(std::shared_ptr<TATensor<double,4>> I993, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I993), ta1_(t2), e0_(e) { }
};

class Task678 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task678(std::shared_ptr<TATensor<double,4>> I993, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I993), ta1_(f1), ta2_(t2) { }
};

class Task679 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task679(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I833), ta1_(t2), e0_(e) { }
};

class Task680 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task680(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I833), ta1_(t2), e0_(e) { }
};

class Task681 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task681(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task682 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task682(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task683 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task683(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task684 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task684(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task685 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task685(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task686 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task686(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task687 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task687(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task688 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task688(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task689 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task689(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task690 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task690(std::shared_ptr<TATensor<double,4>> I833, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I833), ta1_(v2), ta2_(t2) { }
};

class Task691 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task691(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma285, std::shared_ptr<TATensor<double,6>> I841)
   : ta0_(I698), ta1_(Gamma285), ta2_(I841) { }
};

class Task692 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task692(std::shared_ptr<TATensor<double,6>> I841, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I842)
   : ta0_(I841), ta1_(t2), ta2_(I842) { }
};

class Task693 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task693(std::shared_ptr<TATensor<double,4>> I842, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I842), ta1_(f1), ta2_(t2) { }
};

class Task694 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task694(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,3>> Gamma286, std::shared_ptr<TATensor<double,2>> I845)
   : ta0_(I698), ta1_(Gamma286), ta2_(I845) { }
};

class Task695 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task695(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I846)
   : ta0_(I845), ta1_(t2), ta2_(I846) { }
};

class Task696 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task696(std::shared_ptr<TATensor<double,2>> I846, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I846), ta1_(t2), ta2_(f1) { }
};

class Task697 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task697(std::shared_ptr<TATensor<double,2>> I846, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I846), ta1_(t2), ta2_(f1) { }
};

class Task698 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task698(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I900)
   : ta0_(I845), ta1_(t2), ta2_(I900) { }
};

class Task699 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task699(std::shared_ptr<TATensor<double,2>> I900, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I900), ta1_(t2), ta2_(f1) { }
};


}
}
}
#endif
#endif

