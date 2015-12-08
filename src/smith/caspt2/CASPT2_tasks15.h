//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks15.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS15_H
#define __SRC_SMITH_CASPT2_TASKS15_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task700 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task700(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task701 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task701(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task702 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task702(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task703 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task703(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task704 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task704(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task705 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task705(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task706 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task706(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task707 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task707(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma307, std::shared_ptr<TATensor<double,6>> I911)
   : ta0_(I768), ta1_(Gamma307), ta2_(I911) { }
};

class Task708 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task708(std::shared_ptr<TATensor<double,6>> I911, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I912)
   : ta0_(I911), ta1_(t2), ta2_(I912) { }
};

class Task709 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task709(std::shared_ptr<TATensor<double,4>> I912, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I912), ta1_(f1), ta2_(t2) { }
};

class Task710 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task710(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,3>> Gamma308, std::shared_ptr<TATensor<double,2>> I915)
   : ta0_(I768), ta1_(Gamma308), ta2_(I915) { }
};

class Task711 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task711(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I916)
   : ta0_(I915), ta1_(t2), ta2_(I916) { }
};

class Task712 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task712(std::shared_ptr<TATensor<double,2>> I916, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I916), ta1_(t2), ta2_(f1) { }
};

class Task713 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task713(std::shared_ptr<TATensor<double,2>> I916, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I916), ta1_(t2), ta2_(f1) { }
};

class Task714 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task714(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I970)
   : ta0_(I915), ta1_(t2), ta2_(I970) { }
};

class Task715 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task715(std::shared_ptr<TATensor<double,2>> I970, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I970), ta1_(t2), ta2_(f1) { }
};

class Task716 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task716(std::shared_ptr<TATensor<double,2>> I970, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I970), ta1_(t2), ta2_(f1) { }
};

class Task717 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task717(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1021)
   : ta0_(I915), ta1_(t2), ta2_(I1021) { }
};

class Task718 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task718(std::shared_ptr<TATensor<double,2>> I1021, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1021), ta1_(f1), ta2_(t2) { }
};

class Task719 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task719(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1025)
   : ta0_(I915), ta1_(t2), ta2_(I1025) { }
};

class Task720 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task720(std::shared_ptr<TATensor<double,2>> I1025, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1025), ta1_(f1), ta2_(t2) { }
};

class Task721 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task721(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1029)
   : ta0_(I915), ta1_(t2), ta2_(I1029) { }
};

class Task722 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task722(std::shared_ptr<TATensor<double,2>> I1029, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1029), ta1_(f1), ta2_(t2) { }
};

class Task723 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task723(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1033)
   : ta0_(I915), ta1_(t2), ta2_(I1033) { }
};

class Task724 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task724(std::shared_ptr<TATensor<double,2>> I1033, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1033), ta1_(f1), ta2_(t2) { }
};

class Task725 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task725(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I1043)
   : ta0_(I915), ta1_(f1), ta2_(I1043) { }
};

class Task726 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task726(std::shared_ptr<TATensor<double,2>> I1043, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1043), ta1_(t2) { }
};

class Task727 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task727(std::shared_ptr<TATensor<double,2>> I1043, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1043), ta1_(t2) { }
};

class Task728 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task728(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I1075)
   : ta0_(I915), ta1_(f1), ta2_(I1075) { }
};

class Task729 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task729(std::shared_ptr<TATensor<double,2>> I1075, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1075), ta1_(t2) { }
};

class Task730 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task730(std::shared_ptr<TATensor<double,2>> I1075, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1075), ta1_(t2) { }
};

class Task731 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task731(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1089)
   : ta0_(I915), ta1_(t2), ta2_(I1089) { }
};

class Task732 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task732(std::shared_ptr<TATensor<double,4>> I1089, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1089), ta1_(f1), ta2_(t2) { }
};

class Task733 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task733(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1093)
   : ta0_(I915), ta1_(t2), ta2_(I1093) { }
};

class Task734 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task734(std::shared_ptr<TATensor<double,4>> I1093, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1093), ta1_(f1), ta2_(t2) { }
};

class Task735 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task735(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1097)
   : ta0_(I915), ta1_(t2), ta2_(I1097) { }
};

class Task736 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task736(std::shared_ptr<TATensor<double,4>> I1097, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1097), ta1_(f1), ta2_(t2) { }
};

class Task737 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task737(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1101)
   : ta0_(I915), ta1_(t2), ta2_(I1101) { }
};

class Task738 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task738(std::shared_ptr<TATensor<double,4>> I1101, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1101), ta1_(f1), ta2_(t2) { }
};

class Task739 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task739(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1105)
   : ta0_(I915), ta1_(t2), ta2_(I1105) { }
};

class Task740 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task740(std::shared_ptr<TATensor<double,4>> I1105, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1105), ta1_(f1), ta2_(t2) { }
};

class Task741 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task741(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1109)
   : ta0_(I915), ta1_(t2), ta2_(I1109) { }
};

class Task742 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task742(std::shared_ptr<TATensor<double,4>> I1109, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1109), ta1_(f1), ta2_(t2) { }
};

class Task743 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task743(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I915), ta1_(t2), e0_(e) { }
};

class Task744 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task744(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I915), ta1_(t2), e0_(e) { }
};

class Task745 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task745(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I915), ta1_(v2), ta2_(t2) { }
};

class Task746 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task746(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I915), ta1_(v2), ta2_(t2) { }
};

class Task747 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task747(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I915), ta1_(v2), ta2_(t2) { }
};

class Task748 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task748(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I915), ta1_(v2), ta2_(t2) { }
};

class Task749 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task749(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I1279)
   : ta0_(I915), ta1_(h1), ta2_(I1279) { }
};


}
}
}
#endif
#endif

