//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks8.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS8_H
#define __SRC_SMITH_CASPT2_TASKS8_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task350 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task350(std::shared_ptr<TATensor<double,2>> I407, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I408)
   : ta0_(I407), ta1_(t2), ta2_(I408) { }
};

class Task351 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task351(std::shared_ptr<TATensor<double,4>> I408, std::shared_ptr<TATensor<double,4>> Gamma3, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I408), ta1_(Gamma3), ta2_(t2) { }
};

class Task352 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task352(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I410)
   : ta0_(proj), ta1_(I410) { }
};

class Task353 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task353(std::shared_ptr<TATensor<double,2>> I410, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I411)
   : ta0_(I410), ta1_(t2), ta2_(I411) { }
};

class Task354 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task354(std::shared_ptr<TATensor<double,2>> I411, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I411), ta1_(Gamma12), ta2_(t2) { }
};

class Task355 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task355(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I413)
   : ta0_(proj), ta1_(I413) { }
};

class Task356 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task356(std::shared_ptr<TATensor<double,2>> I413, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I414)
   : ta0_(I413), ta1_(t2), ta2_(I414) { }
};

class Task357 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task357(std::shared_ptr<TATensor<double,2>> I414, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I414), ta1_(Gamma12), ta2_(t2) { }
};

class Task358 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task358(std::shared_ptr<TATensor<double,2>> I413, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I570)
   : ta0_(I413), ta1_(t2), ta2_(I570) { }
};

class Task359 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task359(std::shared_ptr<TATensor<double,2>> I570, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I571)
   : ta0_(I570), ta1_(Gamma38), ta2_(I571) { }
};

class Task360 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task360(std::shared_ptr<TATensor<double,4>> I571, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I571), ta1_(t2) { }
};

class Task361 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task361(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I416)
   : ta0_(proj), ta1_(I416) { }
};

class Task362 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task362(std::shared_ptr<TATensor<double,2>> I416, std::shared_ptr<TATensor<double,4>> Gamma152, std::shared_ptr<TATensor<double,2>> I417)
   : ta0_(I416), ta1_(Gamma152), ta2_(I417) { }
};

class Task363 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task363(std::shared_ptr<TATensor<double,2>> I417, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I417), ta1_(t2) { }
};

class Task364 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task364(std::shared_ptr<TATensor<double,2>> I417, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I417), ta1_(t2) { }
};

class Task365 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task365(std::shared_ptr<TATensor<double,2>> I416, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,2>> I626)
   : ta0_(I416), ta1_(Gamma60), ta2_(I626) { }
};

class Task366 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task366(std::shared_ptr<TATensor<double,2>> I626, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I626), ta1_(t2) { }
};

class Task367 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task367(std::shared_ptr<TATensor<double,2>> I626, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I626), ta1_(t2) { }
};

class Task368 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task368(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I422)
   : ta0_(proj), ta1_(I422) { }
};

class Task369 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task369(std::shared_ptr<TATensor<double,2>> I422, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I423)
   : ta0_(I422), ta1_(t2), ta2_(I423) { }
};

class Task370 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task370(std::shared_ptr<TATensor<double,4>> I423, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I423), ta1_(Gamma16), ta2_(t2) { }
};

class Task371 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task371(std::shared_ptr<TATensor<double,2>> I422, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I426)
   : ta0_(I422), ta1_(t2), ta2_(I426) { }
};

class Task372 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task372(std::shared_ptr<TATensor<double,4>> I426, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I426), ta1_(Gamma16), ta2_(t2) { }
};

class Task373 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task373(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I428)
   : ta0_(proj), ta1_(I428) { }
};

class Task374 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task374(std::shared_ptr<TATensor<double,2>> I428, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I429)
   : ta0_(I428), ta1_(t2), ta2_(I429) { }
};

class Task375 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task375(std::shared_ptr<TATensor<double,4>> I429, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I429), ta1_(Gamma16), ta2_(t2) { }
};

class Task376 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task376(std::shared_ptr<TATensor<double,2>> I428, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I435)
   : ta0_(I428), ta1_(t2), ta2_(I435) { }
};

class Task377 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task377(std::shared_ptr<TATensor<double,4>> I435, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I435), ta1_(Gamma16), ta2_(t2) { }
};

class Task378 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task378(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I431)
   : ta0_(proj), ta1_(I431) { }
};

class Task379 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task379(std::shared_ptr<TATensor<double,2>> I431, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I432)
   : ta0_(I431), ta1_(t2), ta2_(I432) { }
};

class Task380 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task380(std::shared_ptr<TATensor<double,4>> I432, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I432), ta1_(Gamma16), ta2_(t2) { }
};

class Task381 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task381(std::shared_ptr<TATensor<double,2>> I431, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I438)
   : ta0_(I431), ta1_(t2), ta2_(I438) { }
};

class Task382 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task382(std::shared_ptr<TATensor<double,4>> I438, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I438), ta1_(Gamma16), ta2_(t2) { }
};

class Task383 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task383(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I440)
   : ta0_(proj), ta1_(I440) { }
};

class Task384 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task384(std::shared_ptr<TATensor<double,2>> I440, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I441)
   : ta0_(I440), ta1_(t2), ta2_(I441) { }
};

class Task385 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task385(std::shared_ptr<TATensor<double,4>> I441, std::shared_ptr<TATensor<double,4>> Gamma22, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I441), ta1_(Gamma22), ta2_(t2) { }
};

class Task386 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task386(std::shared_ptr<TATensor<double,4>> I441, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I441), ta1_(Gamma12), ta2_(t2) { }
};

class Task387 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task387(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I443)
   : ta0_(proj), ta1_(I443) { }
};

class Task388 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task388(std::shared_ptr<TATensor<double,2>> I443, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I444)
   : ta0_(I443), ta1_(t2), ta2_(I444) { }
};

class Task389 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task389(std::shared_ptr<TATensor<double,4>> I444, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> I445)
   : ta0_(I444), ta1_(Gamma12), ta2_(I445) { }
};

class Task390 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task390(std::shared_ptr<TATensor<double,4>> I445, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I445), ta1_(t2) { }
};

class Task391 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task391(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I452)
   : ta0_(proj), ta1_(I452) { }
};

class Task392 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task392(std::shared_ptr<TATensor<double,2>> I452, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I453)
   : ta0_(I452), ta1_(Gamma16), ta2_(I453) { }
};

class Task393 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task393(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I453), ta1_(t2) { }
};

class Task394 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task394(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I453), ta1_(t2) { }
};

class Task395 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task395(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I458)
   : ta0_(proj), ta1_(I458) { }
};

class Task396 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task396(std::shared_ptr<TATensor<double,2>> I458, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I459)
   : ta0_(I458), ta1_(t2), ta2_(I459) { }
};

class Task397 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task397(std::shared_ptr<TATensor<double,4>> I459, std::shared_ptr<TATensor<double,6>> Gamma28, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I459), ta1_(Gamma28), ta2_(t2) { }
};

class Task398 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task398(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I461)
   : ta0_(proj), ta1_(I461) { }
};

class Task399 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task399(std::shared_ptr<TATensor<double,2>> I461, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I462)
   : ta0_(I461), ta1_(t2), ta2_(I462) { }
};


}
}
}
#endif
#endif

