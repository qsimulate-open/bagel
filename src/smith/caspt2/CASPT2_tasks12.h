//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks12.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS12_H
#define __SRC_SMITH_CASPT2_TASKS12_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task550 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task550(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma251, std::shared_ptr<TATensor<double,4>> I710)
   : ta0_(I698), ta1_(Gamma251), ta2_(I710) { }
};

class Task551 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task551(std::shared_ptr<TATensor<double,4>> I710, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I711)
   : ta0_(I710), ta1_(t2), ta2_(I711) { }
};

class Task552 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task552(std::shared_ptr<TATensor<double,4>> I711, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I711), ta1_(t2), ta2_(f1) { }
};

class Task553 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task553(std::shared_ptr<TATensor<double,4>> I710, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I742)
   : ta0_(I710), ta1_(t2), ta2_(I742) { }
};

class Task554 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task554(std::shared_ptr<TATensor<double,4>> I742, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I742), ta1_(f1), ta2_(t2) { }
};

class Task555 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task555(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma252, std::shared_ptr<TATensor<double,6>> I714)
   : ta0_(I698), ta1_(Gamma252), ta2_(I714) { }
};

class Task556 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task556(std::shared_ptr<TATensor<double,6>> I714, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I715)
   : ta0_(I714), ta1_(t2), ta2_(I715) { }
};

class Task557 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task557(std::shared_ptr<TATensor<double,4>> I715, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I715), ta1_(t2), ta2_(f1) { }
};

class Task558 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task558(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma253, std::shared_ptr<TATensor<double,6>> I718)
   : ta0_(I698), ta1_(Gamma253), ta2_(I718) { }
};

class Task559 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task559(std::shared_ptr<TATensor<double,6>> I718, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I718), ta1_(t2) { }
};

class Task560 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task560(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma254, std::shared_ptr<TATensor<double,6>> I721)
   : ta0_(I698), ta1_(Gamma254), ta2_(I721) { }
};

class Task561 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task561(std::shared_ptr<TATensor<double,6>> I721, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I722)
   : ta0_(I721), ta1_(t2), ta2_(I722) { }
};

class Task562 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task562(std::shared_ptr<TATensor<double,4>> I722, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I722), ta1_(f1), ta2_(t2) { }
};

class Task563 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task563(std::shared_ptr<TATensor<double,6>> I721, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I738)
   : ta0_(I721), ta1_(t2), ta2_(I738) { }
};

class Task564 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task564(std::shared_ptr<TATensor<double,4>> I738, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I738), ta1_(t2), ta2_(f1) { }
};

class Task565 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task565(std::shared_ptr<TATensor<double,6>> I721, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I862)
   : ta0_(I721), ta1_(t2), ta2_(I862) { }
};

class Task566 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task566(std::shared_ptr<TATensor<double,4>> I862, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I862), ta1_(t2), e0_(e) { }
};

class Task567 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task567(std::shared_ptr<TATensor<double,4>> I862, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I862), ta1_(f1), ta2_(t2) { }
};

class Task568 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task568(std::shared_ptr<TATensor<double,6>> I721, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I721), ta1_(v2), ta2_(t2) { }
};

class Task569 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task569(std::shared_ptr<TATensor<double,6>> I721, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I721), ta1_(v2), ta2_(t2) { }
};

class Task570 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task570(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma255, std::shared_ptr<TATensor<double,4>> I725)
   : ta0_(I698), ta1_(Gamma255), ta2_(I725) { }
};

class Task571 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task571(std::shared_ptr<TATensor<double,4>> I725, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I726)
   : ta0_(I725), ta1_(t2), ta2_(I726) { }
};

class Task572 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task572(std::shared_ptr<TATensor<double,2>> I726, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I726), ta1_(t2), ta2_(f1) { }
};

class Task573 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task573(std::shared_ptr<TATensor<double,2>> I726, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I726), ta1_(t2), ta2_(f1) { }
};

class Task574 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task574(std::shared_ptr<TATensor<double,4>> I725, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I816)
   : ta0_(I725), ta1_(t2), ta2_(I816) { }
};

class Task575 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task575(std::shared_ptr<TATensor<double,4>> I816, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I816), ta1_(t2), ta2_(f1) { }
};

class Task576 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task576(std::shared_ptr<TATensor<double,4>> I725, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I866)
   : ta0_(I725), ta1_(t2), ta2_(I866) { }
};

class Task577 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task577(std::shared_ptr<TATensor<double,4>> I866, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I866), ta1_(t2), ta2_(f1) { }
};

class Task578 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task578(std::shared_ptr<TATensor<double,4>> I866, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I866), ta1_(t2), ta2_(f1) { }
};

class Task579 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task579(std::shared_ptr<TATensor<double,4>> I725, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I725), ta1_(v2), ta2_(t2) { }
};

class Task580 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task580(std::shared_ptr<TATensor<double,4>> I725, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I725), ta1_(h1), ta2_(t2) { }
};

class Task581 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task581(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma257, std::shared_ptr<TATensor<double,6>> I733)
   : ta0_(I698), ta1_(Gamma257), ta2_(I733) { }
};

class Task582 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task582(std::shared_ptr<TATensor<double,6>> I733, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I734)
   : ta0_(I733), ta1_(t2), ta2_(I734) { }
};

class Task583 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task583(std::shared_ptr<TATensor<double,4>> I734, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I734), ta1_(t2), ta2_(f1) { }
};

class Task584 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task584(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma260, std::shared_ptr<TATensor<double,4>> I745)
   : ta0_(I698), ta1_(Gamma260), ta2_(I745) { }
};

class Task585 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task585(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I746)
   : ta0_(I745), ta1_(t2), ta2_(I746) { }
};

class Task586 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task586(std::shared_ptr<TATensor<double,2>> I746, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I746), ta1_(f1), ta2_(t2) { }
};

class Task587 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task587(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I750)
   : ta0_(I745), ta1_(t2), ta2_(I750) { }
};

class Task588 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task588(std::shared_ptr<TATensor<double,2>> I750, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I750), ta1_(f1), ta2_(t2) { }
};

class Task589 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task589(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I788)
   : ta0_(I745), ta1_(t2), ta2_(I788) { }
};

class Task590 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task590(std::shared_ptr<TATensor<double,4>> I788, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I788), ta1_(f1), ta2_(t2) { }
};

class Task591 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task591(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I792)
   : ta0_(I745), ta1_(t2), ta2_(I792) { }
};

class Task592 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task592(std::shared_ptr<TATensor<double,4>> I792, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I792), ta1_(f1), ta2_(t2) { }
};

class Task593 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task593(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I796)
   : ta0_(I745), ta1_(t2), ta2_(I796) { }
};

class Task594 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task594(std::shared_ptr<TATensor<double,4>> I796, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I796), ta1_(f1), ta2_(t2) { }
};

class Task595 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task595(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I745), ta1_(v2), ta2_(t2) { }
};

class Task596 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task596(std::shared_ptr<TATensor<double,4>> I745, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I745), ta1_(h1), ta2_(t2) { }
};

class Task597 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task597(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,3>> Gamma262, std::shared_ptr<TATensor<double,2>> I753)
   : ta0_(I698), ta1_(Gamma262), ta2_(I753) { }
};

class Task598 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task598(std::shared_ptr<TATensor<double,2>> I753, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I753), ta1_(t2) { }
};

class Task599 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task599(std::shared_ptr<TATensor<double,2>> I753, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I753), ta1_(t2) { }
};


}
}
}
#endif
#endif

