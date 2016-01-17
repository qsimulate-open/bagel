//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks15.h
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
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task700(std::shared_ptr<TATensor<double,2>> I900, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I900), ta1_(t2), ta2_(f1) { }
};

class Task701 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task701(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I951)
   : ta0_(I845), ta1_(t2), ta2_(I951) { }
};

class Task702 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task702(std::shared_ptr<TATensor<double,2>> I951, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I951), ta1_(f1), ta2_(t2) { }
};

class Task703 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task703(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I955)
   : ta0_(I845), ta1_(t2), ta2_(I955) { }
};

class Task704 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task704(std::shared_ptr<TATensor<double,2>> I955, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I955), ta1_(f1), ta2_(t2) { }
};

class Task705 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task705(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I959)
   : ta0_(I845), ta1_(t2), ta2_(I959) { }
};

class Task706 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task706(std::shared_ptr<TATensor<double,2>> I959, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I959), ta1_(f1), ta2_(t2) { }
};

class Task707 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task707(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I963)
   : ta0_(I845), ta1_(t2), ta2_(I963) { }
};

class Task708 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task708(std::shared_ptr<TATensor<double,2>> I963, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I963), ta1_(f1), ta2_(t2) { }
};

class Task709 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task709(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I973)
   : ta0_(I845), ta1_(f1), ta2_(I973) { }
};

class Task710 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task710(std::shared_ptr<TATensor<double,2>> I973, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I973), ta1_(t2) { }
};

class Task711 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task711(std::shared_ptr<TATensor<double,2>> I973, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I973), ta1_(t2) { }
};

class Task712 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task712(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I1005)
   : ta0_(I845), ta1_(f1), ta2_(I1005) { }
};

class Task713 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task713(std::shared_ptr<TATensor<double,2>> I1005, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1005), ta1_(t2) { }
};

class Task714 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task714(std::shared_ptr<TATensor<double,2>> I1005, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1005), ta1_(t2) { }
};

class Task715 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task715(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1019)
   : ta0_(I845), ta1_(t2), ta2_(I1019) { }
};

class Task716 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task716(std::shared_ptr<TATensor<double,4>> I1019, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1019), ta1_(f1), ta2_(t2) { }
};

class Task717 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task717(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1023)
   : ta0_(I845), ta1_(t2), ta2_(I1023) { }
};

class Task718 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task718(std::shared_ptr<TATensor<double,4>> I1023, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1023), ta1_(f1), ta2_(t2) { }
};

class Task719 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task719(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1027)
   : ta0_(I845), ta1_(t2), ta2_(I1027) { }
};

class Task720 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task720(std::shared_ptr<TATensor<double,4>> I1027, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1027), ta1_(f1), ta2_(t2) { }
};

class Task721 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task721(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1031)
   : ta0_(I845), ta1_(t2), ta2_(I1031) { }
};

class Task722 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task722(std::shared_ptr<TATensor<double,4>> I1031, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1031), ta1_(f1), ta2_(t2) { }
};

class Task723 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task723(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1035)
   : ta0_(I845), ta1_(t2), ta2_(I1035) { }
};

class Task724 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task724(std::shared_ptr<TATensor<double,4>> I1035, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1035), ta1_(f1), ta2_(t2) { }
};

class Task725 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task725(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1039)
   : ta0_(I845), ta1_(t2), ta2_(I1039) { }
};

class Task726 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task726(std::shared_ptr<TATensor<double,4>> I1039, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1039), ta1_(f1), ta2_(t2) { }
};

class Task727 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task727(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I845), ta1_(t2), e0_(e) { }
};

class Task728 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task728(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I845), ta1_(t2), e0_(e) { }
};

class Task729 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task729(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I845), ta1_(v2), ta2_(t2) { }
};

class Task730 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task730(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I845), ta1_(v2), ta2_(t2) { }
};

class Task731 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task731(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I845), ta1_(v2), ta2_(t2) { }
};

class Task732 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task732(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I845), ta1_(v2), ta2_(t2) { }
};

class Task733 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task733(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I1209)
   : ta0_(I845), ta1_(h1), ta2_(I1209) { }
};

class Task734 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task734(std::shared_ptr<TATensor<double,4>> I1209, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1209), ta1_(t2) { }
};

class Task735 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task735(std::shared_ptr<TATensor<double,2>> I845, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I1212)
   : ta0_(I845), ta1_(h1), ta2_(I1212) { }
};

class Task736 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task736(std::shared_ptr<TATensor<double,4>> I1212, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1212), ta1_(t2) { }
};

class Task737 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task737(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma299, std::shared_ptr<TATensor<double,6>> I895)
   : ta0_(I698), ta1_(Gamma299), ta2_(I895) { }
};

class Task738 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task738(std::shared_ptr<TATensor<double,6>> I895, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I896)
   : ta0_(I895), ta1_(t2), ta2_(I896) { }
};

class Task739 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task739(std::shared_ptr<TATensor<double,4>> I896, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I896), ta1_(f1), ta2_(t2) { }
};

class Task740 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task740(std::shared_ptr<TATensor<double,6>> I895, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I895), ta1_(v2), ta2_(t2) { }
};

class Task741 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task741(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma304, std::shared_ptr<TATensor<double,6>> I915)
   : ta0_(I698), ta1_(Gamma304), ta2_(I915) { }
};

class Task742 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task742(std::shared_ptr<TATensor<double,6>> I915, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I916)
   : ta0_(I915), ta1_(t2), ta2_(I916) { }
};

class Task743 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task743(std::shared_ptr<TATensor<double,4>> I916, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I916), ta1_(t2), ta2_(f1) { }
};

class Task744 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task744(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma305, std::shared_ptr<TATensor<double,6>> I919)
   : ta0_(I698), ta1_(Gamma305), ta2_(I919) { }
};

class Task745 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task745(std::shared_ptr<TATensor<double,6>> I919, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I920)
   : ta0_(I919), ta1_(t2), ta2_(I920) { }
};

class Task746 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task746(std::shared_ptr<TATensor<double,4>> I920, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I920), ta1_(t2), ta2_(f1) { }
};

class Task747 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task747(std::shared_ptr<TATensor<double,6>> I919, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I919), ta1_(v2), ta2_(t2) { }
};

class Task748 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task748(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma306, std::shared_ptr<TATensor<double,6>> I923)
   : ta0_(I698), ta1_(Gamma306), ta2_(I923) { }
};

class Task749 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task749(std::shared_ptr<TATensor<double,6>> I923, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I923), ta1_(t2) { }
};


}
}
}
#endif
#endif

