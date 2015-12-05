//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks3.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS3_H
#define __SRC_SMITH_RelMRCI_TASKS3_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task100 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task100(std::shared_ptr<TATensor<std::complex<double>,6>> I198, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma63, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I198), ta1_(Gamma63), ta2_(v2) { }
};

class Task101 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task101(std::shared_ptr<TATensor<std::complex<double>,6>> I198, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma64, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I198), ta1_(Gamma64), ta2_(v2) { }
};

class Task102 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task102(std::shared_ptr<TATensor<std::complex<double>,4>> I0, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I204)
   : ta0_(I0), ta1_(v2), ta2_(I204) { }
};

class Task103 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task103(std::shared_ptr<TATensor<std::complex<double>,4>> I204, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma65, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I204), ta1_(Gamma65), ta2_(t2) { }
};

class Task104 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task104(std::shared_ptr<TATensor<std::complex<double>,4>> I0, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I207)
   : ta0_(I0), ta1_(t2), ta2_(I207) { }
};

class Task105 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task105(std::shared_ptr<TATensor<std::complex<double>,4>> I207, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma66, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I207), ta1_(Gamma66), ta2_(v2) { }
};

class Task106 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task106(std::shared_ptr<TATensor<std::complex<double>,4>> I207, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma67, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I207), ta1_(Gamma67), ta2_(v2) { }
};

class Task107 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task107(std::shared_ptr<TATensor<std::complex<double>,4>> I0, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> I213)
   : ta0_(I0), ta1_(Gamma0), ta2_(I213) { }
};

class Task108 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task108(std::shared_ptr<TATensor<std::complex<double>,4>> I213, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I214)
   : ta0_(I213), ta1_(t2), ta2_(I214) { }
};

class Task109 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task109(std::shared_ptr<TATensor<std::complex<double>,4>> I214, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I214), ta1_(v2) { }
};

class Task110 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task110(std::shared_ptr<TATensor<std::complex<double>,4>> I213, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I213), ta1_(t2), ta2_(v2) { }
};

class Task111 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task111(std::shared_ptr<TATensor<std::complex<double>,4>> I0, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma65, std::shared_ptr<TATensor<std::complex<double>,6>> I225)
   : ta0_(I0), ta1_(Gamma65), ta2_(I225) { }
};

class Task112 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task112(std::shared_ptr<TATensor<std::complex<double>,6>> I225, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I225), ta1_(t2), ta2_(v2) { }
};

class Task113 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task113(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I9)
   : ta0_(proj), ta1_(I9) { }
};

class Task114 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task114(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma3, std::shared_ptr<TATensor<std::complex<double>,4>> I10)
   : ta0_(I9), ta1_(Gamma3), ta2_(I10) { }
};

class Task115 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task115(std::shared_ptr<TATensor<std::complex<double>,4>> I10, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I10), ta1_(t2), ta2_(h1) { }
};

class Task116 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task116(std::shared_ptr<TATensor<std::complex<double>,4>> I10, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I10), ta1_(t2), ta2_(v2) { }
};

class Task117 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task117(std::shared_ptr<TATensor<std::complex<double>,4>> I10, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I10), ta1_(t2), ta2_(v2) { }
};

class Task118 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task118(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I13)
   : ta0_(I9), ta1_(h1), ta2_(I13) { }
};

class Task119 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task119(std::shared_ptr<TATensor<std::complex<double>,4>> I13, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I13), ta1_(Gamma4), ta2_(t2) { }
};

class Task120 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task120(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,2>> I16)
   : ta0_(I9), ta1_(Gamma5), ta2_(I16) { }
};

class Task121 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task121(std::shared_ptr<TATensor<std::complex<double>,2>> I16, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I16), ta1_(t2), ta2_(h1) { }
};

class Task122 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task122(std::shared_ptr<TATensor<std::complex<double>,2>> I16, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I16), ta1_(t2), ta2_(h1) { }
};

class Task123 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task123(std::shared_ptr<TATensor<std::complex<double>,2>> I16, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task124 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task124(std::shared_ptr<TATensor<std::complex<double>,2>> I16, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task125 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task125(std::shared_ptr<TATensor<std::complex<double>,2>> I16, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task126 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task126(std::shared_ptr<TATensor<std::complex<double>,2>> I16, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task127 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task127(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> I22)
   : ta0_(I9), ta1_(Gamma4), ta2_(I22) { }
};

class Task128 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task128(std::shared_ptr<TATensor<std::complex<double>,4>> I22, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I22), ta1_(t2), ta2_(h1) { }
};

class Task129 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task129(std::shared_ptr<TATensor<std::complex<double>,4>> I22, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I289)
   : ta0_(I22), ta1_(t2), ta2_(I289) { }
};

class Task130 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task130(std::shared_ptr<TATensor<std::complex<double>,4>> I289, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I289), ta1_(v2) { }
};

class Task131 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task131(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma74, std::shared_ptr<TATensor<std::complex<double>,6>> I231)
   : ta0_(I9), ta1_(Gamma74), ta2_(I231) { }
};

class Task132 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task132(std::shared_ptr<TATensor<std::complex<double>,6>> I231, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I231), ta1_(t2), ta2_(v2) { }
};

class Task133 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task133(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma75, std::shared_ptr<TATensor<std::complex<double>,6>> I234)
   : ta0_(I9), ta1_(Gamma75), ta2_(I234) { }
};

class Task134 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task134(std::shared_ptr<TATensor<std::complex<double>,6>> I234, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I234), ta1_(t2), ta2_(v2) { }
};

class Task135 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task135(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma77, std::shared_ptr<TATensor<std::complex<double>,6>> I240)
   : ta0_(I9), ta1_(Gamma77), ta2_(I240) { }
};

class Task136 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task136(std::shared_ptr<TATensor<std::complex<double>,6>> I240, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I241)
   : ta0_(I240), ta1_(t2), ta2_(I241) { }
};

class Task137 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task137(std::shared_ptr<TATensor<std::complex<double>,4>> I241, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I241), ta1_(v2) { }
};

class Task138 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task138(std::shared_ptr<TATensor<std::complex<double>,6>> I240, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I240), ta1_(t2), ta2_(v2) { }
};

class Task139 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task139(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma78, std::shared_ptr<TATensor<std::complex<double>,6>> I243)
   : ta0_(I9), ta1_(Gamma78), ta2_(I243) { }
};

class Task140 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task140(std::shared_ptr<TATensor<std::complex<double>,6>> I243, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task141 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task141(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma79, std::shared_ptr<TATensor<std::complex<double>,6>> I246)
   : ta0_(I9), ta1_(Gamma79), ta2_(I246) { }
};

class Task142 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task142(std::shared_ptr<TATensor<std::complex<double>,6>> I246, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I246), ta1_(t2), ta2_(v2) { }
};

class Task143 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task143(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma81, std::shared_ptr<TATensor<std::complex<double>,4>> I252)
   : ta0_(I9), ta1_(Gamma81), ta2_(I252) { }
};

class Task144 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task144(std::shared_ptr<TATensor<std::complex<double>,4>> I252, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I253)
   : ta0_(I252), ta1_(t2), ta2_(I253) { }
};

class Task145 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task145(std::shared_ptr<TATensor<std::complex<double>,4>> I253, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I253), ta1_(v2) { }
};

class Task146 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task146(std::shared_ptr<TATensor<std::complex<double>,4>> I252, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I256)
   : ta0_(I252), ta1_(t2), ta2_(I256) { }
};

class Task147 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task147(std::shared_ptr<TATensor<std::complex<double>,4>> I256, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I256), ta1_(v2) { }
};

class Task148 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task148(std::shared_ptr<TATensor<std::complex<double>,4>> I252, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I252), ta1_(t2), ta2_(v2) { }
};

class Task149 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task149(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma84, std::shared_ptr<TATensor<std::complex<double>,4>> I261)
   : ta0_(I9), ta1_(Gamma84), ta2_(I261) { }
};


}
}
}
#endif
#endif

