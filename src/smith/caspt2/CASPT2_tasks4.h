//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks4.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS4_H
#define __SRC_SMITH_CASPT2_TASKS4_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task150 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task150(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I121)
   : ta0_(I120), ta1_(f1), ta2_(I121) { }
};

class Task151 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task151(std::shared_ptr<TATensor<double,4>> I121, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I121), ta1_(Gamma6), ta2_(t2) { }
};

class Task152 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task152(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> I124)
   : ta0_(I120), ta1_(Gamma7), ta2_(I124) { }
};

class Task153 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task153(std::shared_ptr<TATensor<double,4>> I124, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I124), ta1_(v2) { }
};

class Task154 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task154(std::shared_ptr<TATensor<double,4>> I124, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I124), ta1_(t2), ta2_(f1) { }
};

class Task155 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task155(std::shared_ptr<TATensor<double,4>> I124, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I124), ta1_(t2), ta2_(f1) { }
};

class Task156 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task156(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,4>> Gamma34, std::shared_ptr<TATensor<double,4>> I130)
   : ta0_(I120), ta1_(Gamma34), ta2_(I130) { }
};

class Task157 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task157(std::shared_ptr<TATensor<double,4>> I130, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I130), ta1_(t2) { }
};

class Task158 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task158(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I132)
   : ta0_(I120), ta1_(Gamma35), ta2_(I132) { }
};

class Task159 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task159(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2, const double e)
   : ta0_(I132), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task160 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task160(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task161 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task161(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task162 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task162(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task163 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task163(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task164 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task164(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task165 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task165(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task166 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task166(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I146)
   : ta0_(I120), ta1_(f1), ta2_(I146) { }
};

class Task167 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task167(std::shared_ptr<TATensor<double,4>> I146, std::shared_ptr<TATensor<double,6>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I146), ta1_(Gamma51), ta2_(t2) { }
};

class Task168 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task168(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I149)
   : ta0_(I120), ta1_(Gamma38), ta2_(I149) { }
};

class Task169 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task169(std::shared_ptr<TATensor<double,2>> I149, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I149), ta1_(h1) { }
};

class Task170 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task170(std::shared_ptr<TATensor<double,2>> I149, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I149), ta1_(t2), ta2_(f1) { }
};

class Task171 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task171(std::shared_ptr<TATensor<double,2>> I149, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I149), ta1_(t2), ta2_(f1) { }
};

class Task172 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task172(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I160)
   : ta0_(proj), ta1_(I160) { }
};

class Task173 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task173(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma56, std::shared_ptr<TATensor<double,4>> I161)
   : ta0_(I160), ta1_(Gamma56), ta2_(I161) { }
};

class Task174 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task174(std::shared_ptr<TATensor<double,4>> I161, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I161), ta1_(t2), ta2_(f1) { }
};

class Task175 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task175(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma57, std::shared_ptr<TATensor<double,4>> I164)
   : ta0_(I160), ta1_(Gamma57), ta2_(I164) { }
};

class Task176 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task176(std::shared_ptr<TATensor<double,4>> I164, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I164), ta1_(v2) { }
};

class Task177 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task177(std::shared_ptr<TATensor<double,4>> I164, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I164), ta1_(t2), ta2_(f1) { }
};

class Task178 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task178(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma58, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I160), ta1_(Gamma58), ta2_(t2) { }
};

class Task179 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task179(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> I169)
   : ta0_(I160), ta1_(Gamma59), ta2_(I169) { }
};

class Task180 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task180(std::shared_ptr<TATensor<double,4>> I169, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2, const double e)
   : ta0_(I169), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task181 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task181(std::shared_ptr<TATensor<double,4>> I169, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I169), ta1_(t2), ta2_(f1) { }
};

class Task182 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task182(std::shared_ptr<TATensor<double,4>> I169, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I169), ta1_(t2), ta2_(f1) { }
};

class Task183 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task183(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,2>> I172)
   : ta0_(I160), ta1_(Gamma60), ta2_(I172) { }
};

class Task184 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task184(std::shared_ptr<TATensor<double,2>> I172, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I172), ta1_(h1) { }
};

class Task185 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task185(std::shared_ptr<TATensor<double,2>> I172, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I172), ta1_(t2), ta2_(f1) { }
};

class Task186 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task186(std::shared_ptr<TATensor<double,2>> I172, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I172), ta1_(t2), ta2_(f1) { }
};

class Task187 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task187(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I180)
   : ta0_(proj), ta1_(I180) { }
};

class Task188 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task188(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I181)
   : ta0_(I180), ta1_(t2), ta2_(I181) { }
};

class Task189 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task189(std::shared_ptr<TATensor<double,2>> I181, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I181), ta1_(Gamma16), ta2_(f1) { }
};

class Task190 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task190(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I184)
   : ta0_(I180), ta1_(t2), ta2_(I184) { }
};

class Task191 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task191(std::shared_ptr<TATensor<double,2>> I184, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I184), ta1_(Gamma16), ta2_(f1) { }
};

class Task192 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task192(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I187)
   : ta0_(I180), ta1_(f1), ta2_(I187) { }
};

class Task193 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task193(std::shared_ptr<TATensor<double,2>> I187, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I188)
   : ta0_(I187), ta1_(Gamma38), ta2_(I188) { }
};

class Task194 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task194(std::shared_ptr<TATensor<double,4>> I188, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I188), ta1_(t2) { }
};

class Task195 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task195(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I190)
   : ta0_(I180), ta1_(f1), ta2_(I190) { }
};

class Task196 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task196(std::shared_ptr<TATensor<double,2>> I190, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I191)
   : ta0_(I190), ta1_(Gamma38), ta2_(I191) { }
};

class Task197 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task197(std::shared_ptr<TATensor<double,4>> I191, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I191), ta1_(t2) { }
};

class Task198 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task198(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I199)
   : ta0_(I180), ta1_(t2), ta2_(I199) { }
};

class Task199 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task199(std::shared_ptr<TATensor<double,2>> I199, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I199), ta1_(Gamma38), ta2_(f1) { }
};


}
}
}
#endif
#endif

