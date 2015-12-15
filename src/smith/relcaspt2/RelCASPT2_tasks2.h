//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_tasks2.h
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

#ifndef __SRC_SMITH_RelCASPT2_TASKS2_H
#define __SRC_SMITH_RelCASPT2_TASKS2_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelCASPT2{

class Task50 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task50(std::shared_ptr<TATensor<std::complex<double>,4>> I11, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma7, std::shared_ptr<TATensor<std::complex<double>,2>> I20)
   : ta0_(I11), ta1_(Gamma7), ta2_(I20) { }
};

class Task51 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task51(std::shared_ptr<TATensor<std::complex<double>,2>> I20, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I20), ta1_(h1) { }
};

class Task52 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task52(std::shared_ptr<TATensor<std::complex<double>,2>> I20, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I20), ta1_(t2), ta2_(f1) { }
};

class Task53 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task53(std::shared_ptr<TATensor<std::complex<double>,2>> I20, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I20), ta1_(t2), ta2_(f1) { }
};

class Task54 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task54(std::shared_ptr<TATensor<std::complex<double>,4>> I11, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> I26)
   : ta0_(I11), ta1_(Gamma9), ta2_(I26) { }
};

class Task55 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task55(std::shared_ptr<TATensor<std::complex<double>,4>> I26, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I26), ta1_(t2), ta2_(f1) { }
};

class Task56 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task56(std::shared_ptr<TATensor<std::complex<double>,4>> I11, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I11), ta1_(Gamma105), ta2_(v2) { }
};

class Task57 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task57(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I31)
   : ta0_(proj), ta1_(I31) { }
};

class Task58 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task58(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I32)
   : ta0_(I31), ta1_(f1), ta2_(I32) { }
};

class Task59 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task59(std::shared_ptr<TATensor<std::complex<double>,4>> I32, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma3, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I32), ta1_(Gamma3), ta2_(t2) { }
};

class Task60 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task60(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,2>> I35)
   : ta0_(I31), ta1_(f1), ta2_(I35) { }
};

class Task61 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task61(std::shared_ptr<TATensor<std::complex<double>,2>> I35, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma12, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I35), ta1_(Gamma12), ta2_(t2) { }
};

class Task62 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task62(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,2>> I38)
   : ta0_(I31), ta1_(f1), ta2_(I38) { }
};

class Task63 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task63(std::shared_ptr<TATensor<std::complex<double>,2>> I38, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma12, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I38), ta1_(Gamma12), ta2_(t2) { }
};

class Task64 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task64(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma14, std::shared_ptr<TATensor<std::complex<double>,4>> I41)
   : ta0_(I31), ta1_(Gamma14), ta2_(I41) { }
};

class Task65 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task65(std::shared_ptr<TATensor<std::complex<double>,4>> I41, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I41), ta1_(t2) { }
};

class Task66 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task66(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma16, std::shared_ptr<TATensor<std::complex<double>,4>> I45)
   : ta0_(I31), ta1_(Gamma16), ta2_(I45) { }
};

class Task67 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task67(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I45), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task68 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task68(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task69 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task69(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task70 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task70(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task71 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task71(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task72 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task72(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task73 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task73(std::shared_ptr<TATensor<std::complex<double>,4>> I45, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I45), ta1_(t2), ta2_(f1) { }
};

class Task74 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task74(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I63)
   : ta0_(I31), ta1_(f1), ta2_(I63) { }
};

class Task75 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task75(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma22, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I63), ta1_(Gamma22), ta2_(t2) { }
};

class Task76 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task76(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma12, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I63), ta1_(Gamma12), ta2_(t2) { }
};

class Task77 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task77(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I66)
   : ta0_(I31), ta1_(f1), ta2_(I66) { }
};

class Task78 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task78(std::shared_ptr<TATensor<std::complex<double>,4>> I66, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma12, std::shared_ptr<TATensor<std::complex<double>,4>> I67)
   : ta0_(I66), ta1_(Gamma12), ta2_(I67) { }
};

class Task79 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task79(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I67), ta1_(t2) { }
};

class Task80 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task80(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I75)
   : ta0_(I31), ta1_(t2), ta2_(I75) { }
};

class Task81 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task81(std::shared_ptr<TATensor<std::complex<double>,2>> I75, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma16, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I75), ta1_(Gamma16), ta2_(f1) { }
};

class Task82 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task82(std::shared_ptr<TATensor<std::complex<double>,4>> I31, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I78)
   : ta0_(I31), ta1_(t2), ta2_(I78) { }
};

class Task83 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task83(std::shared_ptr<TATensor<std::complex<double>,2>> I78, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma16, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I78), ta1_(Gamma16), ta2_(f1) { }
};

class Task84 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task84(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I80)
   : ta0_(proj), ta1_(I80) { }
};

class Task85 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task85(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I81)
   : ta0_(I80), ta1_(f1), ta2_(I81) { }
};

class Task86 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task86(std::shared_ptr<TATensor<std::complex<double>,4>> I81, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma28, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I81), ta1_(Gamma28), ta2_(t2) { }
};

class Task87 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task87(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma29, std::shared_ptr<TATensor<std::complex<double>,4>> I84)
   : ta0_(I80), ta1_(Gamma29), ta2_(I84) { }
};

class Task88 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task88(std::shared_ptr<TATensor<std::complex<double>,4>> I84, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I84), ta1_(v2) { }
};

class Task89 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task89(std::shared_ptr<TATensor<std::complex<double>,4>> I84, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I84), ta1_(t2), ta2_(f1) { }
};

class Task90 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task90(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma7, std::shared_ptr<TATensor<std::complex<double>,4>> I87)
   : ta0_(I80), ta1_(Gamma7), ta2_(I87) { }
};

class Task91 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task91(std::shared_ptr<TATensor<std::complex<double>,4>> I87, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I87), ta1_(t2), ta2_(f1) { }
};

class Task92 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task92(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma31, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I80), ta1_(Gamma31), ta2_(t2) { }
};

class Task93 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task93(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> I92)
   : ta0_(I80), ta1_(Gamma32), ta2_(I92) { }
};

class Task94 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task94(std::shared_ptr<TATensor<std::complex<double>,4>> I92, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I92), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task95 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task95(std::shared_ptr<TATensor<std::complex<double>,4>> I92, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I92), ta1_(t2), ta2_(f1) { }
};

class Task96 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task96(std::shared_ptr<TATensor<std::complex<double>,4>> I92, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I92), ta1_(t2), ta2_(f1) { }
};

class Task97 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task97(std::shared_ptr<TATensor<std::complex<double>,4>> I92, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I92), ta1_(t2), ta2_(f1) { }
};

class Task98 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task98(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma34, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I80), ta1_(Gamma34), ta2_(t2) { }
};

class Task99 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task99(std::shared_ptr<TATensor<std::complex<double>,4>> I80, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma35, std::shared_ptr<TATensor<std::complex<double>,4>> I100)
   : ta0_(I80), ta1_(Gamma35), ta2_(I100) { }
};


}
}
}
#endif
#endif

