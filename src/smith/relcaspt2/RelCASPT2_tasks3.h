//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_tasks3.h
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

#ifndef __SRC_SMITH_RelCASPT2_TASKS3_H
#define __SRC_SMITH_RelCASPT2_TASKS3_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelCASPT2{

class Task100 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task100(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma35, std::shared_ptr<TATensor<std::complex<double>,4>> I100)
   : ta0_(I80), ta1_(Gamma35), ta2_(I100) { }
};

class Task101 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task101(std::shared_ptr<TATensor<std::complex<double>,4>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I100), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task102 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task102(std::shared_ptr<TATensor<std::complex<double>,4>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I100), ta1_(t2), ta2_(f1) { }
};

class Task103 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task103(std::shared_ptr<TATensor<std::complex<double>,4>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I100), ta1_(t2), ta2_(f1) { }
};

class Task104 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task104(std::shared_ptr<TATensor<std::complex<double>,4>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I100), ta1_(t2), ta2_(f1) { }
};

class Task105 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task105(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I106)
   : ta0_(I80), ta1_(f1), ta2_(I106) { }
};

class Task106 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task106(std::shared_ptr<TATensor<std::complex<double>,4>> I106, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma37, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I106), ta1_(Gamma37), ta2_(t2) { }
};

class Task107 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task107(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,2>> I109)
   : ta0_(I80), ta1_(Gamma38), ta2_(I109) { }
};

class Task108 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task108(std::shared_ptr<TATensor<std::complex<double>,2>> I109, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I109), ta1_(h1) { }
};

class Task109 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task109(std::shared_ptr<TATensor<std::complex<double>,2>> I109, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I109), ta1_(t2), ta2_(f1) { }
};

class Task110 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task110(std::shared_ptr<TATensor<std::complex<double>,2>> I109, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I109), ta1_(t2), ta2_(f1) { }
};

class Task111 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task111(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I120)
   : ta0_(proj), ta1_(I120) { }
};

class Task112 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task112(std::shared_ptr<TATensor<std::complex<double>,4>> I120, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I121)
   : ta0_(I120), ta1_(f1), ta2_(I121) { }
};

class Task113 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task113(std::shared_ptr<TATensor<std::complex<double>,4>> I121, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma6, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I121), ta1_(Gamma6), ta2_(t2) { }
};

class Task114 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task114(std::shared_ptr<TATensor<std::complex<double>,4>> I120, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma7, std::shared_ptr<TATensor<std::complex<double>,4>> I124)
   : ta0_(I120), ta1_(Gamma7), ta2_(I124) { }
};

class Task115 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task115(std::shared_ptr<TATensor<std::complex<double>,4>> I124, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I124), ta1_(v2) { }
};

class Task116 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task116(std::shared_ptr<TATensor<std::complex<double>,4>> I124, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I124), ta1_(t2), ta2_(f1) { }
};

class Task117 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task117(std::shared_ptr<TATensor<std::complex<double>,4>> I124, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I124), ta1_(t2), ta2_(f1) { }
};

class Task118 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task118(std::shared_ptr<TATensor<std::complex<double>,4>> I120, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma34, std::shared_ptr<TATensor<std::complex<double>,4>> I130)
   : ta0_(I120), ta1_(Gamma34), ta2_(I130) { }
};

class Task119 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task119(std::shared_ptr<TATensor<std::complex<double>,4>> I130, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I130), ta1_(t2) { }
};

class Task120 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task120(std::shared_ptr<TATensor<std::complex<double>,4>> I120, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma35, std::shared_ptr<TATensor<std::complex<double>,4>> I132)
   : ta0_(I120), ta1_(Gamma35), ta2_(I132) { }
};

class Task121 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task121(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I132), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task122 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task122(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task123 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task123(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task124 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task124(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task125 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task125(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task126 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task126(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task127 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task127(std::shared_ptr<TATensor<std::complex<double>,4>> I132, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task128 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task128(std::shared_ptr<TATensor<std::complex<double>,4>> I120, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I146)
   : ta0_(I120), ta1_(f1), ta2_(I146) { }
};

class Task129 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task129(std::shared_ptr<TATensor<std::complex<double>,4>> I146, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma51, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I146), ta1_(Gamma51), ta2_(t2) { }
};

class Task130 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task130(std::shared_ptr<TATensor<std::complex<double>,4>> I120, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,2>> I149)
   : ta0_(I120), ta1_(Gamma38), ta2_(I149) { }
};

class Task131 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task131(std::shared_ptr<TATensor<std::complex<double>,2>> I149, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I149), ta1_(h1) { }
};

class Task132 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task132(std::shared_ptr<TATensor<std::complex<double>,2>> I149, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I149), ta1_(t2), ta2_(f1) { }
};

class Task133 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task133(std::shared_ptr<TATensor<std::complex<double>,2>> I149, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I149), ta1_(t2), ta2_(f1) { }
};

class Task134 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task134(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I160)
   : ta0_(proj), ta1_(I160) { }
};

class Task135 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task135(std::shared_ptr<TATensor<std::complex<double>,4>> I160, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma56, std::shared_ptr<TATensor<std::complex<double>,4>> I161)
   : ta0_(I160), ta1_(Gamma56), ta2_(I161) { }
};

class Task136 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task136(std::shared_ptr<TATensor<std::complex<double>,4>> I161, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I161), ta1_(t2), ta2_(f1) { }
};

class Task137 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task137(std::shared_ptr<TATensor<std::complex<double>,4>> I160, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma57, std::shared_ptr<TATensor<std::complex<double>,4>> I164)
   : ta0_(I160), ta1_(Gamma57), ta2_(I164) { }
};

class Task138 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task138(std::shared_ptr<TATensor<std::complex<double>,4>> I164, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I164), ta1_(v2) { }
};

class Task139 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task139(std::shared_ptr<TATensor<std::complex<double>,4>> I164, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I164), ta1_(t2), ta2_(f1) { }
};

class Task140 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task140(std::shared_ptr<TATensor<std::complex<double>,4>> I160, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma58, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I160), ta1_(Gamma58), ta2_(t2) { }
};

class Task141 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task141(std::shared_ptr<TATensor<std::complex<double>,4>> I160, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma59, std::shared_ptr<TATensor<std::complex<double>,4>> I169)
   : ta0_(I160), ta1_(Gamma59), ta2_(I169) { }
};

class Task142 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task142(std::shared_ptr<TATensor<std::complex<double>,4>> I169, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I169), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task143 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task143(std::shared_ptr<TATensor<std::complex<double>,4>> I169, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I169), ta1_(t2), ta2_(f1) { }
};

class Task144 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task144(std::shared_ptr<TATensor<std::complex<double>,4>> I169, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I169), ta1_(t2), ta2_(f1) { }
};

class Task145 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task145(std::shared_ptr<TATensor<std::complex<double>,4>> I160, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,2>> I172)
   : ta0_(I160), ta1_(Gamma60), ta2_(I172) { }
};

class Task146 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task146(std::shared_ptr<TATensor<std::complex<double>,2>> I172, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I172), ta1_(h1) { }
};

class Task147 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task147(std::shared_ptr<TATensor<std::complex<double>,2>> I172, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I172), ta1_(t2), ta2_(f1) { }
};

class Task148 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task148(std::shared_ptr<TATensor<std::complex<double>,2>> I172, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I172), ta1_(t2), ta2_(f1) { }
};

class Task149 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task149(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I180)
   : ta0_(proj), ta1_(I180) { }
};


}
}
}
#endif
#endif

