//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks3.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS3_H
#define __SRC_SMITH_CASPT2_TASKS3_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task100 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task100(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> Gamma14, std::shared_ptr<TATensor<double,4>> I41)
   : ta0_(I31), ta1_(Gamma14), ta2_(I41) { }
};

class Task101 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task101(std::shared_ptr<TATensor<double,4>> I41, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I41), ta1_(t2) { }
};

class Task102 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task102(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> I45)
   : ta0_(I31), ta1_(Gamma16), ta2_(I45) { }
};

class Task103 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task103(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I45), ta1_(t2), e0_(e) { }
};

class Task104 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task104(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task105 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task105(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task106 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task106(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task107 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task107(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task108 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task108(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task109 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task109(std::shared_ptr<TATensor<double,4>> I45, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task110 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task110(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I63)
   : ta0_(I31), ta1_(f1), ta2_(I63) { }
};

class Task111 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task111(std::shared_ptr<TATensor<double,4>> I63, std::shared_ptr<TATensor<double,4>> Gamma22, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I63), ta1_(Gamma22), ta2_(t2) { }
};

class Task112 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task112(std::shared_ptr<TATensor<double,4>> I63, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I63), ta1_(Gamma12), ta2_(t2) { }
};

class Task113 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task113(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I66)
   : ta0_(I31), ta1_(f1), ta2_(I66) { }
};

class Task114 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task114(std::shared_ptr<TATensor<double,4>> I66, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> I67)
   : ta0_(I66), ta1_(Gamma12), ta2_(I67) { }
};

class Task115 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task115(std::shared_ptr<TATensor<double,4>> I67, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I67), ta1_(t2) { }
};

class Task116 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task116(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I75)
   : ta0_(I31), ta1_(t2), ta2_(I75) { }
};

class Task117 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task117(std::shared_ptr<TATensor<double,2>> I75, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I75), ta1_(Gamma16), ta2_(f1) { }
};

class Task118 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task118(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I78)
   : ta0_(I31), ta1_(t2), ta2_(I78) { }
};

class Task119 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task119(std::shared_ptr<TATensor<double,2>> I78, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I78), ta1_(Gamma16), ta2_(f1) { }
};

class Task120 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task120(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I80)
   : ta0_(proj), ta1_(I80) { }
};

class Task121 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task121(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I81)
   : ta0_(I80), ta1_(f1), ta2_(I81) { }
};

class Task122 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task122(std::shared_ptr<TATensor<double,4>> I81, std::shared_ptr<TATensor<double,6>> Gamma28, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I81), ta1_(Gamma28), ta2_(t2) { }
};

class Task123 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task123(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I84)
   : ta0_(I80), ta1_(Gamma29), ta2_(I84) { }
};

class Task124 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task124(std::shared_ptr<TATensor<double,4>> I84, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I84), ta1_(t2), ta2_(f1) { }
};

class Task125 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task125(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> I87)
   : ta0_(I80), ta1_(Gamma7), ta2_(I87) { }
};

class Task126 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task126(std::shared_ptr<TATensor<double,4>> I87, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I87), ta1_(t2), ta2_(f1) { }
};

class Task127 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task127(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,4>> Gamma31, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I80), ta1_(Gamma31), ta2_(t2) { }
};

class Task128 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task128(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> I92)
   : ta0_(I80), ta1_(Gamma32), ta2_(I92) { }
};

class Task129 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task129(std::shared_ptr<TATensor<double,4>> I92, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I92), ta1_(t2), e0_(e) { }
};

class Task130 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task130(std::shared_ptr<TATensor<double,4>> I92, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I92), ta1_(t2), ta2_(f1) { }
};

class Task131 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task131(std::shared_ptr<TATensor<double,4>> I92, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I92), ta1_(t2), ta2_(f1) { }
};

class Task132 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task132(std::shared_ptr<TATensor<double,4>> I92, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I92), ta1_(t2), ta2_(f1) { }
};

class Task133 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task133(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,4>> Gamma34, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I80), ta1_(Gamma34), ta2_(t2) { }
};

class Task134 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task134(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I100)
   : ta0_(I80), ta1_(Gamma35), ta2_(I100) { }
};

class Task135 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task135(std::shared_ptr<TATensor<double,4>> I100, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I100), ta1_(t2), e0_(e) { }
};

class Task136 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task136(std::shared_ptr<TATensor<double,4>> I100, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I100), ta1_(t2), ta2_(f1) { }
};

class Task137 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task137(std::shared_ptr<TATensor<double,4>> I100, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I100), ta1_(t2), ta2_(f1) { }
};

class Task138 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task138(std::shared_ptr<TATensor<double,4>> I100, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I100), ta1_(t2), ta2_(f1) { }
};

class Task139 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task139(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I106)
   : ta0_(I80), ta1_(f1), ta2_(I106) { }
};

class Task140 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task140(std::shared_ptr<TATensor<double,4>> I106, std::shared_ptr<TATensor<double,6>> Gamma37, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I106), ta1_(Gamma37), ta2_(t2) { }
};

class Task141 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task141(std::shared_ptr<TATensor<double,4>> I80, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I109)
   : ta0_(I80), ta1_(Gamma38), ta2_(I109) { }
};

class Task142 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task142(std::shared_ptr<TATensor<double,2>> I109, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I109), ta1_(t2), ta2_(f1) { }
};

class Task143 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task143(std::shared_ptr<TATensor<double,2>> I109, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I109), ta1_(t2), ta2_(f1) { }
};

class Task144 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task144(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I120)
   : ta0_(proj), ta1_(I120) { }
};

class Task145 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task145(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I121)
   : ta0_(I120), ta1_(f1), ta2_(I121) { }
};

class Task146 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task146(std::shared_ptr<TATensor<double,4>> I121, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I121), ta1_(Gamma6), ta2_(t2) { }
};

class Task147 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task147(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> I124)
   : ta0_(I120), ta1_(Gamma7), ta2_(I124) { }
};

class Task148 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task148(std::shared_ptr<TATensor<double,4>> I124, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I124), ta1_(t2), ta2_(f1) { }
};

class Task149 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task149(std::shared_ptr<TATensor<double,4>> I124, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I124), ta1_(t2), ta2_(f1) { }
};


}
}
}
#endif
#endif

