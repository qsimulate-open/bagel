//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks4.h
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

#ifndef __SRC_SMITH_MRCI_TASKS4_H
#define __SRC_SMITH_MRCI_TASKS4_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task150 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task150(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,2>> I16)
   : ta0_(I9), ta1_(Gamma5), ta2_(I16) { }
};

class Task151 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task151(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I16), ta1_(t2), ta2_(h1) { }
};

class Task152 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task152(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I16), ta1_(t2), ta2_(h1) { }
};

class Task153 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task153(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task154 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task154(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task155 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task155(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task156 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task156(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};

class Task157 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task157(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma7, std::shared_ptr<TATensor<double,4>> I22)
   : ta0_(I9), ta1_(Gamma7), ta2_(I22) { }
};

class Task158 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task158(std::shared_ptr<TATensor<double,4>> I22, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I22), ta1_(t2), ta2_(h1) { }
};

class Task159 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task159(std::shared_ptr<TATensor<double,4>> I22, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I22), ta1_(t2), ta2_(v2) { }
};

class Task160 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task160(std::shared_ptr<TATensor<double,4>> I22, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I22), ta1_(t2), ta2_(v2) { }
};

class Task161 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task161(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma97, std::shared_ptr<TATensor<double,6>> I300)
   : ta0_(I9), ta1_(Gamma97), ta2_(I300) { }
};

class Task162 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task162(std::shared_ptr<TATensor<double,6>> I300, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I300), ta1_(t2), ta2_(v2) { }
};

class Task163 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task163(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma98, std::shared_ptr<TATensor<double,6>> I303)
   : ta0_(I9), ta1_(Gamma98), ta2_(I303) { }
};

class Task164 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task164(std::shared_ptr<TATensor<double,6>> I303, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I303), ta1_(t2), ta2_(v2) { }
};

class Task165 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task165(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma100, std::shared_ptr<TATensor<double,6>> I309)
   : ta0_(I9), ta1_(Gamma100), ta2_(I309) { }
};

class Task166 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task166(std::shared_ptr<TATensor<double,6>> I309, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I310)
   : ta0_(I309), ta1_(t2), ta2_(I310) { }
};

class Task167 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task167(std::shared_ptr<TATensor<double,4>> I310, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I310), ta1_(v2) { }
};

class Task168 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task168(std::shared_ptr<TATensor<double,6>> I309, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I309), ta1_(t2), ta2_(v2) { }
};

class Task169 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task169(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma101, std::shared_ptr<TATensor<double,6>> I312)
   : ta0_(I9), ta1_(Gamma101), ta2_(I312) { }
};

class Task170 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task170(std::shared_ptr<TATensor<double,6>> I312, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I312), ta1_(t2), ta2_(v2) { }
};

class Task171 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task171(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma102, std::shared_ptr<TATensor<double,6>> I315)
   : ta0_(I9), ta1_(Gamma102), ta2_(I315) { }
};

class Task172 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task172(std::shared_ptr<TATensor<double,6>> I315, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I315), ta1_(t2), ta2_(v2) { }
};

class Task173 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task173(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma104, std::shared_ptr<TATensor<double,4>> I321)
   : ta0_(I9), ta1_(Gamma104), ta2_(I321) { }
};

class Task174 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task174(std::shared_ptr<TATensor<double,4>> I321, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I322)
   : ta0_(I321), ta1_(t2), ta2_(I322) { }
};

class Task175 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task175(std::shared_ptr<TATensor<double,4>> I322, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I322), ta1_(v2) { }
};

class Task176 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task176(std::shared_ptr<TATensor<double,4>> I321, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I325)
   : ta0_(I321), ta1_(t2), ta2_(I325) { }
};

class Task177 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task177(std::shared_ptr<TATensor<double,4>> I325, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I325), ta1_(v2) { }
};

class Task178 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task178(std::shared_ptr<TATensor<double,4>> I321, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I321), ta1_(t2), ta2_(v2) { }
};

class Task179 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task179(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma107, std::shared_ptr<TATensor<double,4>> I330)
   : ta0_(I9), ta1_(Gamma107), ta2_(I330) { }
};

class Task180 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task180(std::shared_ptr<TATensor<double,4>> I330, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I330), ta1_(t2), ta2_(v2) { }
};

class Task181 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task181(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma109, std::shared_ptr<TATensor<double,4>> I336)
   : ta0_(I9), ta1_(Gamma109), ta2_(I336) { }
};

class Task182 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task182(std::shared_ptr<TATensor<double,4>> I336, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I336), ta1_(t2), ta2_(v2) { }
};

class Task183 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task183(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma114, std::shared_ptr<TATensor<double,6>> I351)
   : ta0_(I9), ta1_(Gamma114), ta2_(I351) { }
};

class Task184 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task184(std::shared_ptr<TATensor<double,6>> I351, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I351), ta1_(t2), ta2_(v2) { }
};

class Task185 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task185(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma115, std::shared_ptr<TATensor<double,6>> I354)
   : ta0_(I9), ta1_(Gamma115), ta2_(I354) { }
};

class Task186 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task186(std::shared_ptr<TATensor<double,6>> I354, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I354), ta1_(t2), ta2_(v2) { }
};

class Task187 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task187(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma119, std::shared_ptr<TATensor<double,6>> I366)
   : ta0_(I9), ta1_(Gamma119), ta2_(I366) { }
};

class Task188 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task188(std::shared_ptr<TATensor<double,6>> I366, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I366), ta1_(t2), ta2_(v2) { }
};

class Task189 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task189(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,8>> Gamma122, std::shared_ptr<TATensor<double,6>> I375)
   : ta0_(I9), ta1_(Gamma122), ta2_(I375) { }
};

class Task190 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task190(std::shared_ptr<TATensor<double,6>> I375, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I375), ta1_(t2), ta2_(v2) { }
};

class Task191 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task191(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma547, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I9), ta1_(Gamma547), ta2_(t2) { }
};

class Task192 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task192(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma548, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I9), ta1_(Gamma548), ta2_(t2) { }
};

class Task193 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task193(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I27)
   : ta0_(proj), ta1_(I27) { }
};

class Task194 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task194(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I28)
   : ta0_(I27), ta1_(h1), ta2_(I28) { }
};

class Task195 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task195(std::shared_ptr<TATensor<double,4>> I28, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I28), ta1_(Gamma2), ta2_(t2) { }
};

class Task196 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task196(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I31)
   : ta0_(I27), ta1_(h1), ta2_(I31) { }
};

class Task197 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task197(std::shared_ptr<TATensor<double,2>> I31, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I31), ta1_(Gamma10), ta2_(t2) { }
};

class Task198 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task198(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I34)
   : ta0_(I27), ta1_(h1), ta2_(I34) { }
};

class Task199 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task199(std::shared_ptr<TATensor<double,2>> I34, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I34), ta1_(Gamma10), ta2_(t2) { }
};


}
}
}
#endif
#endif

