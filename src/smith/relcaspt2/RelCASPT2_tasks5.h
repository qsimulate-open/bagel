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
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task200(std::shared_ptr<TATensor<std::complex<double>,4>> I290, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I290), ta1_(Gamma32), ta2_(v2) { }
};

class Task201 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task201(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I292)
   : ta0_(proj), ta1_(I292) { }
};

class Task202 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task202(std::shared_ptr<TATensor<std::complex<double>,4>> I292, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I292), ta1_(Gamma38), ta2_(h1) { }
};

class Task203 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task203(std::shared_ptr<TATensor<std::complex<double>,4>> I292, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma35, std::shared_ptr<TATensor<std::complex<double>,4>> I315)
   : ta0_(I292), ta1_(Gamma35), ta2_(I315) { }
};

class Task204 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task204(std::shared_ptr<TATensor<std::complex<double>,4>> I315, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I315), ta1_(v2) { }
};

class Task205 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task205(std::shared_ptr<TATensor<std::complex<double>,4>> I292, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma7, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I292), ta1_(Gamma7), ta2_(v2) { }
};

class Task206 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task206(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I294)
   : ta0_(proj), ta1_(I294) { }
};

class Task207 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task207(std::shared_ptr<TATensor<std::complex<double>,4>> I294, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I294), ta1_(Gamma60), ta2_(h1) { }
};

class Task208 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task208(std::shared_ptr<TATensor<std::complex<double>,4>> I294, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma59, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I294), ta1_(Gamma59), ta2_(v2) { }
};

class Task209 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task209(std::shared_ptr<TATensor<std::complex<double>,4>> I294, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma57, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I294), ta1_(Gamma57), ta2_(v2) { }
};

class Task210 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task210(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I296)
   : ta0_(proj), ta1_(I296) { }
};

class Task211 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task211(std::shared_ptr<TATensor<std::complex<double>,4>> I296, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma92, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I296), ta1_(Gamma92), ta2_(v2) { }
};

class Task212 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task212(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I302)
   : ta0_(proj), ta1_(I302) { }
};

class Task213 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task213(std::shared_ptr<TATensor<std::complex<double>,4>> I302, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma16, std::shared_ptr<TATensor<std::complex<double>,4>> I303)
   : ta0_(I302), ta1_(Gamma16), ta2_(I303) { }
};

class Task214 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task214(std::shared_ptr<TATensor<std::complex<double>,4>> I303, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I303), ta1_(v2) { }
};

class Task215 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task215(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I326)
   : ta0_(proj), ta1_(I326) { }
};

class Task216 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task216(std::shared_ptr<TATensor<std::complex<double>,4>> I326, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I326), ta1_(v2) { }
};

class Task217 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task217(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I328)
   : ta0_(proj), ta1_(I328) { }
};

class Task218 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task218(std::shared_ptr<TATensor<std::complex<double>,4>> I328, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma38, std::shared_ptr<TATensor<std::complex<double>,4>> I329)
   : ta0_(I328), ta1_(Gamma38), ta2_(I329) { }
};

class Task219 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task219(std::shared_ptr<TATensor<std::complex<double>,4>> I329, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I329), ta1_(v2) { }
};

class Task220 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task220(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I332)
   : ta0_(proj), ta1_(I332) { }
};

class Task221 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task221(std::shared_ptr<TATensor<std::complex<double>,4>> I332, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma60, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I332), ta1_(Gamma60), ta2_(v2) { }
};


}
}
}
#endif
#endif

