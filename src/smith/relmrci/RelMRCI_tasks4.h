//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks4.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS4_H
#define __SRC_SMITH_RelMRCI_TASKS4_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task150 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task150(std::shared_ptr<TATensor<std::complex<double>,4>> I267, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I267), ta1_(t2), ta2_(v2) { }
};

class Task151 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task151(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma92, std::shared_ptr<TATensor<std::complex<double>,6>> I285)
   : ta0_(I9), ta1_(Gamma92), ta2_(I285) { }
};

class Task152 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task152(std::shared_ptr<TATensor<std::complex<double>,6>> I285, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I285), ta1_(t2), ta2_(v2) { }
};

class Task153 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task153(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma95, std::shared_ptr<TATensor<std::complex<double>,6>> I294)
   : ta0_(I9), ta1_(Gamma95), ta2_(I294) { }
};

class Task154 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task154(std::shared_ptr<TATensor<std::complex<double>,6>> I294, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I294), ta1_(t2), ta2_(v2) { }
};

class Task155 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task155(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma98, std::shared_ptr<TATensor<std::complex<double>,4>> I303)
   : ta0_(I9), ta1_(Gamma98), ta2_(I303) { }
};

class Task156 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task156(std::shared_ptr<TATensor<std::complex<double>,4>> I303, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I303), ta1_(t2), ta2_(v2) { }
};

class Task157 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task157(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma411, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I9), ta1_(Gamma411), ta2_(t2) { }
};

class Task158 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task158(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma412, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I9), ta1_(Gamma412), ta2_(t2) { }
};

class Task159 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task159(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I24)
   : ta0_(proj), ta1_(I24) { }
};

class Task160 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task160(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I25)
   : ta0_(I24), ta1_(h1), ta2_(I25) { }
};

class Task161 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task161(std::shared_ptr<TATensor<std::complex<double>,4>> I25, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I25), ta1_(Gamma2), ta2_(t2) { }
};

class Task162 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task162(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,2>> I28)
   : ta0_(I24), ta1_(h1), ta2_(I28) { }
};

class Task163 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task163(std::shared_ptr<TATensor<std::complex<double>,2>> I28, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I28), ta1_(Gamma9), ta2_(t2) { }
};

class Task164 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task164(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,2>> I31)
   : ta0_(I24), ta1_(h1), ta2_(I31) { }
};

class Task165 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task165(std::shared_ptr<TATensor<std::complex<double>,2>> I31, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I31), ta1_(Gamma9), ta2_(t2) { }
};

class Task166 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task166(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I34)
   : ta0_(I24), ta1_(Gamma11), ta2_(I34) { }
};

class Task167 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task167(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I34), ta1_(t2), ta2_(h1) { }
};

class Task168 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task168(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I34), ta1_(t2), ta2_(h1) { }
};

class Task169 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task169(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I34), ta1_(t2), ta2_(h1) { }
};

class Task170 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task170(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I34), ta1_(t2), ta2_(h1) { }
};

class Task171 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task171(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I34), ta1_(t2), ta2_(h1) { }
};

class Task172 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task172(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I34), ta1_(t2), ta2_(h1) { }
};

class Task173 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task173(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task174 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task174(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task175 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task175(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task176 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task176(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task177 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task177(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task178 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task178(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task179 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task179(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task180 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task180(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task181 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task181(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task182 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task182(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task183 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task183(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I502)
   : ta0_(I34), ta1_(t2), ta2_(I502) { }
};

class Task184 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task184(std::shared_ptr<TATensor<std::complex<double>,4>> I502, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I502), ta1_(v2) { }
};

class Task185 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task185(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I505)
   : ta0_(I34), ta1_(t2), ta2_(I505) { }
};

class Task186 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task186(std::shared_ptr<TATensor<std::complex<double>,4>> I505, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I505), ta1_(v2) { }
};

class Task187 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task187(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I514)
   : ta0_(I34), ta1_(t2), ta2_(I514) { }
};

class Task188 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task188(std::shared_ptr<TATensor<std::complex<double>,4>> I514, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I514), ta1_(v2) { }
};

class Task189 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task189(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I517)
   : ta0_(I34), ta1_(t2), ta2_(I517) { }
};

class Task190 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task190(std::shared_ptr<TATensor<std::complex<double>,4>> I517, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I517), ta1_(v2) { }
};

class Task191 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task191(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task192 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task192(std::shared_ptr<TATensor<std::complex<double>,4>> I34, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I34), ta1_(t2), ta2_(v2) { }
};

class Task193 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task193(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I52)
   : ta0_(I24), ta1_(h1), ta2_(I52) { }
};

class Task194 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task194(std::shared_ptr<TATensor<std::complex<double>,4>> I52, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I52), ta1_(Gamma9), ta2_(t2) { }
};

class Task195 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task195(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I55)
   : ta0_(I24), ta1_(h1), ta2_(I55) { }
};

class Task196 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task196(std::shared_ptr<TATensor<std::complex<double>,4>> I55, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I55), ta1_(Gamma9), ta2_(t2) { }
};

class Task197 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task197(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I58)
   : ta0_(I24), ta1_(t2), ta2_(I58) { }
};

class Task198 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task198(std::shared_ptr<TATensor<std::complex<double>,2>> I58, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I58), ta1_(Gamma11), ta2_(h1) { }
};

class Task199 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task199(std::shared_ptr<TATensor<std::complex<double>,2>> I58, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I58), ta1_(Gamma160), ta2_(v2) { }
};


}
}
}
#endif
#endif

