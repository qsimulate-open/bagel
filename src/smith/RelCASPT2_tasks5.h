//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_tasks5.h
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

#ifndef __SRC_SMITH_RelCASPT2_TASKS5_H
#define __SRC_SMITH_RelCASPT2_TASKS5_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelCASPT2{

class Task200 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task200(std::shared_ptr<TATensor<std::complex<double>,4>> I265, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I266)
   : ta0_(I265), ta1_(f1), ta2_(I266) { }
};

class Task201 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task201(std::shared_ptr<TATensor<std::complex<double>,4>> I266, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma59, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I266), ta1_(Gamma59), ta2_(t2) { }
};

class Task202 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task202(std::shared_ptr<TATensor<std::complex<double>,4>> I265, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,4>> I269)
   : ta0_(I265), ta1_(Gamma60), ta2_(I269) { }
};

class Task203 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task203(std::shared_ptr<TATensor<std::complex<double>,4>> I269, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I269), ta1_(t2), ta2_(f1) { }
};

class Task204 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task204(std::shared_ptr<TATensor<std::complex<double>,4>> I269, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I269), ta1_(t2), ta2_(f1) { }
};

class Task205 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task205(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I271)
   : ta0_(proj), ta1_(I271) { }
};

class Task206 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task206(std::shared_ptr<TATensor<std::complex<double>,4>> I271, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma92, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I271), ta1_(Gamma92), ta2_(t2) { }
};

class Task207 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task207(std::shared_ptr<TATensor<std::complex<double>,4>> I271, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,4>> I301)
   : ta0_(I271), ta1_(Gamma60), ta2_(I301) { }
};

class Task208 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task208(std::shared_ptr<TATensor<std::complex<double>,4>> I301, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I301), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task209 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task209(std::shared_ptr<TATensor<std::complex<double>,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task210 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task210(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I348)
   : ta0_(proj), ta1_(I348) { }
};

class Task211 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task211(std::shared_ptr<TATensor<std::complex<double>,4>> I348, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma7, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I348), ta1_(Gamma7), ta2_(h1) { }
};

class Task212 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task212(std::shared_ptr<TATensor<std::complex<double>,4>> I348, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma107, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I348), ta1_(Gamma107), ta2_(v2) { }
};

class Task213 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task213(std::shared_ptr<TATensor<std::complex<double>,4>> I348, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma6, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I348), ta1_(Gamma6), ta2_(v2) { }
};

class Task214 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task214(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I350)
   : ta0_(proj), ta1_(I350) { }
};

class Task215 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task215(std::shared_ptr<TATensor<std::complex<double>,4>> I350, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I350), ta1_(Gamma38), ta2_(h1) { }
};

class Task216 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task216(std::shared_ptr<TATensor<std::complex<double>,4>> I350, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma35, std::shared_ptr<TATensor<std::complex<double>,4>> I367)
   : ta0_(I350), ta1_(Gamma35), ta2_(I367) { }
};

class Task217 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task217(std::shared_ptr<TATensor<std::complex<double>,4>> I367, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I367), ta1_(v2) { }
};

class Task218 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task218(std::shared_ptr<TATensor<std::complex<double>,4>> I350, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma29, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I350), ta1_(Gamma29), ta2_(v2) { }
};

class Task219 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task219(std::shared_ptr<TATensor<std::complex<double>,4>> I350, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I350), ta1_(Gamma32), ta2_(v2) { }
};

class Task220 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task220(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I352)
   : ta0_(proj), ta1_(I352) { }
};

class Task221 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task221(std::shared_ptr<TATensor<std::complex<double>,4>> I352, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I352), ta1_(Gamma38), ta2_(h1) { }
};

class Task222 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task222(std::shared_ptr<TATensor<std::complex<double>,4>> I352, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma35, std::shared_ptr<TATensor<std::complex<double>,4>> I375)
   : ta0_(I352), ta1_(Gamma35), ta2_(I375) { }
};

class Task223 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task223(std::shared_ptr<TATensor<std::complex<double>,4>> I375, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I375), ta1_(v2) { }
};

class Task224 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task224(std::shared_ptr<TATensor<std::complex<double>,4>> I352, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma7, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I352), ta1_(Gamma7), ta2_(v2) { }
};

class Task225 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task225(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I354)
   : ta0_(proj), ta1_(I354) { }
};

class Task226 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task226(std::shared_ptr<TATensor<std::complex<double>,4>> I354, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I354), ta1_(Gamma60), ta2_(h1) { }
};

class Task227 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task227(std::shared_ptr<TATensor<std::complex<double>,4>> I354, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma59, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I354), ta1_(Gamma59), ta2_(v2) { }
};

class Task228 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task228(std::shared_ptr<TATensor<std::complex<double>,4>> I354, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma57, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I354), ta1_(Gamma57), ta2_(v2) { }
};

class Task229 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task229(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I356)
   : ta0_(proj), ta1_(I356) { }
};

class Task230 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task230(std::shared_ptr<TATensor<std::complex<double>,4>> I356, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma94, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I356), ta1_(Gamma94), ta2_(v2) { }
};

class Task231 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task231(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I362)
   : ta0_(proj), ta1_(I362) { }
};

class Task232 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task232(std::shared_ptr<TATensor<std::complex<double>,4>> I362, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma16, std::shared_ptr<TATensor<std::complex<double>,4>> I363)
   : ta0_(I362), ta1_(Gamma16), ta2_(I363) { }
};

class Task233 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task233(std::shared_ptr<TATensor<std::complex<double>,4>> I363, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I363), ta1_(v2) { }
};

class Task234 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task234(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I386)
   : ta0_(proj), ta1_(I386) { }
};

class Task235 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task235(std::shared_ptr<TATensor<std::complex<double>,4>> I386, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I386), ta1_(v2) { }
};

class Task236 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task236(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I388)
   : ta0_(proj), ta1_(I388) { }
};

class Task237 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task237(std::shared_ptr<TATensor<std::complex<double>,4>> I388, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,4>> I389)
   : ta0_(I388), ta1_(Gamma38), ta2_(I389) { }
};

class Task238 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task238(std::shared_ptr<TATensor<std::complex<double>,4>> I389, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I389), ta1_(v2) { }
};

class Task239 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task239(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I392)
   : ta0_(proj), ta1_(I392) { }
};

class Task240 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task240(std::shared_ptr<TATensor<std::complex<double>,4>> I392, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I392), ta1_(Gamma60), ta2_(v2) { }
};


}
}
}
#endif
#endif

