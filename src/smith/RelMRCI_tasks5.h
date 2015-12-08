//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks5.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS5_H
#define __SRC_SMITH_RelMRCI_TASKS5_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task200 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task200(std::shared_ptr<TATensor<std::complex<double>,2>> I58, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I58), ta1_(Gamma160), ta2_(v2) { }
};

class Task201 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task201(std::shared_ptr<TATensor<std::complex<double>,2>> I58, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I58), ta1_(Gamma9), ta2_(v2) { }
};

class Task202 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task202(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I61)
   : ta0_(I24), ta1_(t2), ta2_(I61) { }
};

class Task203 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task203(std::shared_ptr<TATensor<std::complex<double>,2>> I61, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I61), ta1_(Gamma11), ta2_(h1) { }
};

class Task204 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task204(std::shared_ptr<TATensor<std::complex<double>,2>> I61, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I61), ta1_(Gamma160), ta2_(v2) { }
};

class Task205 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task205(std::shared_ptr<TATensor<std::complex<double>,2>> I61, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I61), ta1_(Gamma9), ta2_(v2) { }
};

class Task206 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task206(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I306)
   : ta0_(I24), ta1_(t2), ta2_(I306) { }
};

class Task207 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task207(std::shared_ptr<TATensor<std::complex<double>,4>> I306, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma99, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I306), ta1_(Gamma99), ta2_(v2) { }
};

class Task208 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task208(std::shared_ptr<TATensor<std::complex<double>,4>> I306, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma66, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I306), ta1_(Gamma66), ta2_(v2) { }
};

class Task209 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task209(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I312)
   : ta0_(I24), ta1_(v2), ta2_(I312) { }
};

class Task210 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task210(std::shared_ptr<TATensor<std::complex<double>,4>> I312, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I312), ta1_(Gamma0), ta2_(t2) { }
};

class Task211 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task211(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I315)
   : ta0_(I24), ta1_(v2), ta2_(I315) { }
};

class Task212 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task212(std::shared_ptr<TATensor<std::complex<double>,4>> I315, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I315), ta1_(Gamma0), ta2_(t2) { }
};

class Task213 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task213(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I318)
   : ta0_(I24), ta1_(v2), ta2_(I318) { }
};

class Task214 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task214(std::shared_ptr<TATensor<std::complex<double>,4>> I318, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I318), ta1_(Gamma0), ta2_(t2) { }
};

class Task215 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task215(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I321)
   : ta0_(I24), ta1_(v2), ta2_(I321) { }
};

class Task216 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task216(std::shared_ptr<TATensor<std::complex<double>,4>> I321, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I321), ta1_(Gamma2), ta2_(t2) { }
};

class Task217 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task217(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I324)
   : ta0_(I24), ta1_(v2), ta2_(I324) { }
};

class Task218 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task218(std::shared_ptr<TATensor<std::complex<double>,4>> I324, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I324), ta1_(Gamma105), ta2_(t2) { }
};

class Task219 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task219(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I327)
   : ta0_(I24), ta1_(v2), ta2_(I327) { }
};

class Task220 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task220(std::shared_ptr<TATensor<std::complex<double>,4>> I327, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I327), ta1_(Gamma105), ta2_(t2) { }
};

class Task221 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task221(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I330)
   : ta0_(I24), ta1_(v2), ta2_(I330) { }
};

class Task222 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task222(std::shared_ptr<TATensor<std::complex<double>,4>> I330, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma1, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I330), ta1_(Gamma1), ta2_(t2) { }
};

class Task223 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task223(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I333)
   : ta0_(I24), ta1_(v2), ta2_(I333) { }
};

class Task224 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task224(std::shared_ptr<TATensor<std::complex<double>,4>> I333, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma65, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I333), ta1_(Gamma65), ta2_(t2) { }
};

class Task225 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task225(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I336)
   : ta0_(I24), ta1_(v2), ta2_(I336) { }
};

class Task226 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task226(std::shared_ptr<TATensor<std::complex<double>,4>> I336, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I336), ta1_(Gamma105), ta2_(t2) { }
};

class Task227 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task227(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I339)
   : ta0_(I24), ta1_(v2), ta2_(I339) { }
};

class Task228 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task228(std::shared_ptr<TATensor<std::complex<double>,4>> I339, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma110, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I339), ta1_(Gamma110), ta2_(t2) { }
};

class Task229 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task229(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I342)
   : ta0_(I24), ta1_(v2), ta2_(I342) { }
};

class Task230 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task230(std::shared_ptr<TATensor<std::complex<double>,4>> I342, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I342), ta1_(Gamma105), ta2_(t2) { }
};

class Task231 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task231(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I345)
   : ta0_(I24), ta1_(v2), ta2_(I345) { }
};

class Task232 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task232(std::shared_ptr<TATensor<std::complex<double>,4>> I345, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I345), ta1_(Gamma105), ta2_(t2) { }
};

class Task233 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task233(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I348)
   : ta0_(I24), ta1_(v2), ta2_(I348) { }
};

class Task234 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task234(std::shared_ptr<TATensor<std::complex<double>,2>> I348, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I348), ta1_(Gamma9), ta2_(t2) { }
};

class Task235 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task235(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I351)
   : ta0_(I24), ta1_(v2), ta2_(I351) { }
};

class Task236 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task236(std::shared_ptr<TATensor<std::complex<double>,2>> I351, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I351), ta1_(Gamma9), ta2_(t2) { }
};

class Task237 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task237(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I354)
   : ta0_(I24), ta1_(t2), ta2_(I354) { }
};

class Task238 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task238(std::shared_ptr<TATensor<std::complex<double>,4>> I354, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> I355)
   : ta0_(I354), ta1_(Gamma160), ta2_(I355) { }
};

class Task239 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task239(std::shared_ptr<TATensor<std::complex<double>,4>> I355, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I355), ta1_(v2) { }
};

class Task240 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task240(std::shared_ptr<TATensor<std::complex<double>,4>> I354, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I354), ta1_(Gamma0), ta2_(v2) { }
};

class Task241 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task241(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I357)
   : ta0_(I24), ta1_(t2), ta2_(I357) { }
};

class Task242 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task242(std::shared_ptr<TATensor<std::complex<double>,4>> I357, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> I358)
   : ta0_(I357), ta1_(Gamma160), ta2_(I358) { }
};

class Task243 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task243(std::shared_ptr<TATensor<std::complex<double>,4>> I358, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I358), ta1_(v2) { }
};

class Task244 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task244(std::shared_ptr<TATensor<std::complex<double>,4>> I357, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I357), ta1_(Gamma2), ta2_(v2) { }
};

class Task245 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task245(std::shared_ptr<TATensor<std::complex<double>,4>> I357, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma128, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I357), ta1_(Gamma128), ta2_(v2) { }
};

class Task246 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task246(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I360)
   : ta0_(I24), ta1_(t2), ta2_(I360) { }
};

class Task247 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task247(std::shared_ptr<TATensor<std::complex<double>,4>> I360, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> I361)
   : ta0_(I360), ta1_(Gamma160), ta2_(I361) { }
};

class Task248 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task248(std::shared_ptr<TATensor<std::complex<double>,4>> I361, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I361), ta1_(v2) { }
};

class Task249 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task249(std::shared_ptr<TATensor<std::complex<double>,4>> I360, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I360), ta1_(Gamma2), ta2_(v2) { }
};


}
}
}
#endif
#endif

