//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks5.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS5_H
#define __SRC_SMITH_CASPT2_TASKS5_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task200 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task200(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I217)
   : ta0_(I204), ta1_(f1), ta2_(I217) { }
};

class Task201 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task201(std::shared_ptr<TATensor<double,2>> I217, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I217), ta1_(Gamma60), ta2_(t2) { }
};

class Task202 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task202(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I220)
   : ta0_(I204), ta1_(f1), ta2_(I220) { }
};

class Task203 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task203(std::shared_ptr<TATensor<double,2>> I220, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I220), ta1_(Gamma60), ta2_(t2) { }
};

class Task204 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task204(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I223)
   : ta0_(I204), ta1_(t2), ta2_(I223) { }
};

class Task205 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task205(std::shared_ptr<TATensor<double,2>> I223, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I223), ta1_(Gamma38), ta2_(f1) { }
};

class Task206 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task206(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I226)
   : ta0_(I204), ta1_(t2), ta2_(I226) { }
};

class Task207 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task207(std::shared_ptr<TATensor<double,2>> I226, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I226), ta1_(Gamma38), ta2_(f1) { }
};

class Task208 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task208(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> Gamma79, std::shared_ptr<TATensor<double,4>> I229)
   : ta0_(I204), ta1_(Gamma79), ta2_(I229) { }
};

class Task209 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task209(std::shared_ptr<TATensor<double,4>> I229, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I229), ta1_(t2) { }
};

class Task210 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task210(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I233)
   : ta0_(I204), ta1_(Gamma38), ta2_(I233) { }
};

class Task211 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task211(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I233), ta1_(t2), e0_(e) { }
};

class Task212 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task212(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I233), ta1_(t2), ta2_(f1) { }
};

class Task213 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task213(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I233), ta1_(t2), ta2_(f1) { }
};

class Task214 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task214(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I233), ta1_(t2), ta2_(f1) { }
};

class Task215 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task215(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I233), ta1_(t2), ta2_(f1) { }
};

class Task216 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task216(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I233), ta1_(t2), ta2_(f1) { }
};

class Task217 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task217(std::shared_ptr<TATensor<double,4>> I233, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I233), ta1_(t2), ta2_(f1) { }
};

class Task218 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task218(std::shared_ptr<TATensor<double,4>> I204, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I251)
   : ta0_(I204), ta1_(f1), ta2_(I251) { }
};

class Task219 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task219(std::shared_ptr<TATensor<double,4>> I251, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I251), ta1_(Gamma60), ta2_(t2) { }
};

class Task220 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task220(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I253)
   : ta0_(proj), ta1_(I253) { }
};

class Task221 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task221(std::shared_ptr<TATensor<double,4>> I253, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I254)
   : ta0_(I253), ta1_(f1), ta2_(I254) { }
};

class Task222 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task222(std::shared_ptr<TATensor<double,4>> I254, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I254), ta1_(Gamma59), ta2_(t2) { }
};

class Task223 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task223(std::shared_ptr<TATensor<double,4>> I253, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> I257)
   : ta0_(I253), ta1_(Gamma60), ta2_(I257) { }
};

class Task224 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task224(std::shared_ptr<TATensor<double,4>> I257, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I257), ta1_(t2), ta2_(f1) { }
};

class Task225 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task225(std::shared_ptr<TATensor<double,4>> I257, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I257), ta1_(t2), ta2_(f1) { }
};

class Task226 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task226(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I259)
   : ta0_(proj), ta1_(I259) { }
};

class Task227 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task227(std::shared_ptr<TATensor<double,4>> I259, std::shared_ptr<TATensor<double,4>> Gamma90, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I259), ta1_(Gamma90), ta2_(t2) { }
};

class Task228 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task228(std::shared_ptr<TATensor<double,4>> I259, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I259), ta1_(Gamma60), ta2_(t2), e0_(e) { }
};

class Task229 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task229(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task230 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task230(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I288)
   : ta0_(proj), ta1_(I288) { }
};

class Task231 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task231(std::shared_ptr<TATensor<double,4>> I288, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I288), ta1_(Gamma92), ta2_(v2) { }
};

class Task232 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task232(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I290)
   : ta0_(proj), ta1_(I290) { }
};

class Task233 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task233(std::shared_ptr<TATensor<double,4>> I290, std::shared_ptr<TATensor<double,6>> Gamma105, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I290), ta1_(Gamma105), ta2_(v2) { }
};

class Task234 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task234(std::shared_ptr<TATensor<double,4>> I290, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I290), ta1_(Gamma6), ta2_(v2) { }
};

class Task235 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task235(std::shared_ptr<TATensor<double,4>> I290, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I290), ta1_(Gamma7), ta2_(h1) { }
};

class Task236 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task236(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I294)
   : ta0_(proj), ta1_(I294) { }
};

class Task237 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task237(std::shared_ptr<TATensor<double,4>> I294, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> I295)
   : ta0_(I294), ta1_(Gamma16), ta2_(I295) { }
};

class Task238 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task238(std::shared_ptr<TATensor<double,4>> I295, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I295), ta1_(v2) { }
};

class Task239 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task239(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I298)
   : ta0_(proj), ta1_(I298) { }
};

class Task240 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task240(std::shared_ptr<TATensor<double,4>> I298, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I299)
   : ta0_(I298), ta1_(Gamma35), ta2_(I299) { }
};

class Task241 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task241(std::shared_ptr<TATensor<double,4>> I299, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I299), ta1_(v2) { }
};

class Task242 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task242(std::shared_ptr<TATensor<double,4>> I298, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I298), ta1_(Gamma29), ta2_(v2) { }
};

class Task243 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task243(std::shared_ptr<TATensor<double,4>> I298, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I298), ta1_(Gamma32), ta2_(v2) { }
};

class Task244 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task244(std::shared_ptr<TATensor<double,4>> I298, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I298), ta1_(Gamma38), ta2_(h1) { }
};

class Task245 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task245(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I306)
   : ta0_(proj), ta1_(I306) { }
};

class Task246 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task246(std::shared_ptr<TATensor<double,4>> I306, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I307)
   : ta0_(I306), ta1_(Gamma35), ta2_(I307) { }
};

class Task247 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task247(std::shared_ptr<TATensor<double,4>> I307, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I307), ta1_(v2) { }
};

class Task248 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task248(std::shared_ptr<TATensor<double,4>> I306, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I306), ta1_(Gamma7), ta2_(v2) { }
};

class Task249 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task249(std::shared_ptr<TATensor<double,4>> I306, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I306), ta1_(Gamma38), ta2_(h1) { }
};


}
}
}
#endif
#endif

