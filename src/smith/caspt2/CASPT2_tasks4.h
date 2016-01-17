//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks4.h
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
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task150(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,4>> Gamma34, std::shared_ptr<TATensor<double,4>> I130)
   : ta0_(I120), ta1_(Gamma34), ta2_(I130) { }
};

class Task151 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task151(std::shared_ptr<TATensor<double,4>> I130, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I130), ta1_(t2) { }
};

class Task152 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task152(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I132)
   : ta0_(I120), ta1_(Gamma35), ta2_(I132) { }
};

class Task153 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task153(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I132), ta1_(t2), e0_(e) { }
};

class Task154 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task154(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task155 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task155(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task156 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task156(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task157 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task157(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task158 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task158(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task159 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task159(std::shared_ptr<TATensor<double,4>> I132, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I132), ta1_(t2), ta2_(f1) { }
};

class Task160 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task160(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I146)
   : ta0_(I120), ta1_(f1), ta2_(I146) { }
};

class Task161 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task161(std::shared_ptr<TATensor<double,4>> I146, std::shared_ptr<TATensor<double,6>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I146), ta1_(Gamma51), ta2_(t2) { }
};

class Task162 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task162(std::shared_ptr<TATensor<double,4>> I120, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I149)
   : ta0_(I120), ta1_(Gamma38), ta2_(I149) { }
};

class Task163 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task163(std::shared_ptr<TATensor<double,2>> I149, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I149), ta1_(t2), ta2_(f1) { }
};

class Task164 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task164(std::shared_ptr<TATensor<double,2>> I149, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I149), ta1_(t2), ta2_(f1) { }
};

class Task165 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task165(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I160)
   : ta0_(proj), ta1_(I160) { }
};

class Task166 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task166(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma56, std::shared_ptr<TATensor<double,4>> I161)
   : ta0_(I160), ta1_(Gamma56), ta2_(I161) { }
};

class Task167 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task167(std::shared_ptr<TATensor<double,4>> I161, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I161), ta1_(t2), ta2_(f1) { }
};

class Task168 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task168(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma57, std::shared_ptr<TATensor<double,4>> I164)
   : ta0_(I160), ta1_(Gamma57), ta2_(I164) { }
};

class Task169 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task169(std::shared_ptr<TATensor<double,4>> I164, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I164), ta1_(t2), ta2_(f1) { }
};

class Task170 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task170(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma58, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I160), ta1_(Gamma58), ta2_(t2) { }
};

class Task171 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task171(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> I169)
   : ta0_(I160), ta1_(Gamma59), ta2_(I169) { }
};

class Task172 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task172(std::shared_ptr<TATensor<double,4>> I169, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I169), ta1_(t2), e0_(e) { }
};

class Task173 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task173(std::shared_ptr<TATensor<double,4>> I169, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I169), ta1_(t2), ta2_(f1) { }
};

class Task174 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task174(std::shared_ptr<TATensor<double,4>> I169, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I169), ta1_(t2), ta2_(f1) { }
};

class Task175 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task175(std::shared_ptr<TATensor<double,4>> I160, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,2>> I172)
   : ta0_(I160), ta1_(Gamma60), ta2_(I172) { }
};

class Task176 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task176(std::shared_ptr<TATensor<double,2>> I172, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I172), ta1_(t2), ta2_(f1) { }
};

class Task177 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task177(std::shared_ptr<TATensor<double,2>> I172, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I172), ta1_(t2), ta2_(f1) { }
};

class Task178 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task178(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I180)
   : ta0_(proj), ta1_(I180) { }
};

class Task179 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task179(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I181)
   : ta0_(I180), ta1_(t2), ta2_(I181) { }
};

class Task180 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task180(std::shared_ptr<TATensor<double,2>> I181, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I181), ta1_(Gamma16), ta2_(f1) { }
};

class Task181 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task181(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I184)
   : ta0_(I180), ta1_(t2), ta2_(I184) { }
};

class Task182 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task182(std::shared_ptr<TATensor<double,2>> I184, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I184), ta1_(Gamma16), ta2_(f1) { }
};

class Task183 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task183(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I187)
   : ta0_(I180), ta1_(f1), ta2_(I187) { }
};

class Task184 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task184(std::shared_ptr<TATensor<double,2>> I187, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I188)
   : ta0_(I187), ta1_(Gamma38), ta2_(I188) { }
};

class Task185 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task185(std::shared_ptr<TATensor<double,4>> I188, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I188), ta1_(t2) { }
};

class Task186 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task186(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I190)
   : ta0_(I180), ta1_(f1), ta2_(I190) { }
};

class Task187 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task187(std::shared_ptr<TATensor<double,2>> I190, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I191)
   : ta0_(I190), ta1_(Gamma38), ta2_(I191) { }
};

class Task188 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task188(std::shared_ptr<TATensor<double,4>> I191, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I191), ta1_(t2) { }
};

class Task189 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task189(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I199)
   : ta0_(I180), ta1_(t2), ta2_(I199) { }
};

class Task190 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task190(std::shared_ptr<TATensor<double,2>> I199, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I199), ta1_(Gamma38), ta2_(f1) { }
};

class Task191 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task191(std::shared_ptr<TATensor<double,4>> I180, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I202)
   : ta0_(I180), ta1_(t2), ta2_(I202) { }
};

class Task192 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task192(std::shared_ptr<TATensor<double,2>> I202, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I202), ta1_(Gamma38), ta2_(f1) { }
};

class Task193 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task193(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I204)
   : ta0_(proj), ta1_(I204) { }
};

class Task194 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task194(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I205)
   : ta0_(I204), ta1_(f1), ta2_(I205) { }
};

class Task195 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task195(std::shared_ptr<TATensor<double,4>> I205, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I206)
   : ta0_(I205), ta1_(Gamma35), ta2_(I206) { }
};

class Task196 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task196(std::shared_ptr<TATensor<double,4>> I206, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I206), ta1_(t2) { }
};

class Task197 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task197(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I208)
   : ta0_(I204), ta1_(f1), ta2_(I208) { }
};

class Task198 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task198(std::shared_ptr<TATensor<double,4>> I208, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I208), ta1_(Gamma32), ta2_(t2) { }
};

class Task199 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task199(std::shared_ptr<TATensor<double,4>> I208, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I208), ta1_(Gamma35), ta2_(t2) { }
};


}
}
}
#endif
#endif

