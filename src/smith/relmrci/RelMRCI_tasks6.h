//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks6.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS6_H
#define __SRC_SMITH_RelMRCI_TASKS6_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task250 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task250(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I363)
   : ta0_(I24), ta1_(t2), ta2_(I363) { }
};

class Task251 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task251(std::shared_ptr<TATensor<std::complex<double>,4>> I363, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> I364)
   : ta0_(I363), ta1_(Gamma160), ta2_(I364) { }
};

class Task252 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task252(std::shared_ptr<TATensor<std::complex<double>,4>> I364, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I364), ta1_(v2) { }
};

class Task253 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task253(std::shared_ptr<TATensor<std::complex<double>,4>> I363, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I363), ta1_(Gamma2), ta2_(v2) { }
};

class Task254 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task254(std::shared_ptr<TATensor<std::complex<double>,4>> I363, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma128, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I363), ta1_(Gamma128), ta2_(v2) { }
};

class Task255 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task255(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I366)
   : ta0_(I24), ta1_(t2), ta2_(I366) { }
};

class Task256 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task256(std::shared_ptr<TATensor<std::complex<double>,4>> I366, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> I367)
   : ta0_(I366), ta1_(Gamma160), ta2_(I367) { }
};

class Task257 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task257(std::shared_ptr<TATensor<std::complex<double>,4>> I367, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I367), ta1_(v2) { }
};

class Task258 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task258(std::shared_ptr<TATensor<std::complex<double>,4>> I366, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I366), ta1_(Gamma2), ta2_(v2) { }
};

class Task259 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task259(std::shared_ptr<TATensor<std::complex<double>,4>> I366, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma128, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I366), ta1_(Gamma128), ta2_(v2) { }
};

class Task260 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task260(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I369)
   : ta0_(I24), ta1_(t2), ta2_(I369) { }
};

class Task261 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task261(std::shared_ptr<TATensor<std::complex<double>,4>> I369, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> I370)
   : ta0_(I369), ta1_(Gamma160), ta2_(I370) { }
};

class Task262 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task262(std::shared_ptr<TATensor<std::complex<double>,4>> I370, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I370), ta1_(v2) { }
};

class Task263 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task263(std::shared_ptr<TATensor<std::complex<double>,4>> I369, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I369), ta1_(Gamma0), ta2_(v2) { }
};

class Task264 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task264(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I456)
   : ta0_(I24), ta1_(t2), ta2_(I456) { }
};

class Task265 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task265(std::shared_ptr<TATensor<std::complex<double>,4>> I456, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I456), ta1_(Gamma105), ta2_(v2) { }
};

class Task266 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task266(std::shared_ptr<TATensor<std::complex<double>,4>> I456, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma151, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I456), ta1_(Gamma151), ta2_(v2) { }
};

class Task267 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task267(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I459)
   : ta0_(I24), ta1_(t2), ta2_(I459) { }
};

class Task268 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task268(std::shared_ptr<TATensor<std::complex<double>,4>> I459, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma105, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I459), ta1_(Gamma105), ta2_(v2) { }
};

class Task269 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task269(std::shared_ptr<TATensor<std::complex<double>,4>> I459, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma151, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I459), ta1_(Gamma151), ta2_(v2) { }
};

class Task270 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task270(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I468)
   : ta0_(I24), ta1_(v2), ta2_(I468) { }
};

class Task271 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task271(std::shared_ptr<TATensor<std::complex<double>,4>> I468, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I468), ta1_(Gamma9), ta2_(t2) { }
};

class Task272 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task272(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I471)
   : ta0_(I24), ta1_(v2), ta2_(I471) { }
};

class Task273 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task273(std::shared_ptr<TATensor<std::complex<double>,4>> I471, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I471), ta1_(Gamma9), ta2_(t2) { }
};

class Task274 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task274(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I474)
   : ta0_(I24), ta1_(v2), ta2_(I474) { }
};

class Task275 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task275(std::shared_ptr<TATensor<std::complex<double>,4>> I474, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I474), ta1_(Gamma9), ta2_(t2) { }
};

class Task276 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task276(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I477)
   : ta0_(I24), ta1_(v2), ta2_(I477) { }
};

class Task277 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task277(std::shared_ptr<TATensor<std::complex<double>,4>> I477, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I477), ta1_(Gamma9), ta2_(t2) { }
};

class Task278 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task278(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I480)
   : ta0_(I24), ta1_(v2), ta2_(I480) { }
};

class Task279 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task279(std::shared_ptr<TATensor<std::complex<double>,4>> I480, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I480), ta1_(Gamma9), ta2_(t2) { }
};

class Task280 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task280(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I483)
   : ta0_(I24), ta1_(v2), ta2_(I483) { }
};

class Task281 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task281(std::shared_ptr<TATensor<std::complex<double>,4>> I483, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I483), ta1_(Gamma9), ta2_(t2) { }
};

class Task282 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task282(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I486)
   : ta0_(I24), ta1_(v2), ta2_(I486) { }
};

class Task283 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task283(std::shared_ptr<TATensor<std::complex<double>,4>> I486, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma159, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I486), ta1_(Gamma159), ta2_(t2) { }
};

class Task284 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task284(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I531)
   : ta0_(I24), ta1_(t2), ta2_(I531) { }
};

class Task285 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task285(std::shared_ptr<TATensor<std::complex<double>,4>> I531, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma174, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I531), ta1_(Gamma174), ta2_(v2) { }
};

class Task286 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task286(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I534)
   : ta0_(I24), ta1_(t2), ta2_(I534) { }
};

class Task287 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task287(std::shared_ptr<TATensor<std::complex<double>,4>> I534, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I534), ta1_(Gamma9), ta2_(v2) { }
};

class Task288 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task288(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I537)
   : ta0_(I24), ta1_(t2), ta2_(I537) { }
};

class Task289 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task289(std::shared_ptr<TATensor<std::complex<double>,4>> I537, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma174, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I537), ta1_(Gamma174), ta2_(v2) { }
};

class Task290 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task290(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I540)
   : ta0_(I24), ta1_(t2), ta2_(I540) { }
};

class Task291 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task291(std::shared_ptr<TATensor<std::complex<double>,4>> I540, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma174, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I540), ta1_(Gamma174), ta2_(v2) { }
};

class Task292 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task292(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma413, std::shared_ptr<TATensor<std::complex<double>,4>> I1264)
   : ta0_(I24), ta1_(Gamma413), ta2_(I1264) { }
};

class Task293 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task293(std::shared_ptr<TATensor<std::complex<double>,4>> I1264, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1264), ta1_(t2) { }
};

class Task294 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task294(std::shared_ptr<TATensor<std::complex<double>,4>> I24, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma415, std::shared_ptr<TATensor<std::complex<double>,4>> I1268)
   : ta0_(I24), ta1_(Gamma415), ta2_(I1268) { }
};

class Task295 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task295(std::shared_ptr<TATensor<std::complex<double>,4>> I1268, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1268), ta1_(t2) { }
};

class Task296 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task296(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I63)
   : ta0_(proj), ta1_(I63) { }
};

class Task297 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task297(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I64)
   : ta0_(I63), ta1_(h1), ta2_(I64) { }
};

class Task298 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task298(std::shared_ptr<TATensor<std::complex<double>,4>> I64, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I64), ta1_(Gamma4), ta2_(t2) { }
};

class Task299 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task299(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,4>> I67)
   : ta0_(I63), ta1_(Gamma5), ta2_(I67) { }
};


}
}
}
#endif
#endif

