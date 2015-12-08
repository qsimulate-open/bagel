//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks12.h
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

#ifndef __SRC_SMITH_MRCI_TASKS12_H
#define __SRC_SMITH_MRCI_TASKS12_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task550 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task550(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> Gamma570, std::shared_ptr<TATensor<double,4>> I1716)
   : ta0_(I108), ta1_(Gamma570), ta2_(I1716) { }
};

class Task551 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task551(std::shared_ptr<TATensor<double,4>> I1716, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1716), ta1_(t2) { }
};

class Task552 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task552(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I144)
   : ta0_(proj), ta1_(I144) { }
};

class Task553 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task553(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma48, std::shared_ptr<TATensor<double,4>> I145)
   : ta0_(I144), ta1_(Gamma48), ta2_(I145) { }
};

class Task554 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task554(std::shared_ptr<TATensor<double,4>> I145, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I145), ta1_(t2), ta2_(h1) { }
};

class Task555 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task555(std::shared_ptr<TATensor<double,4>> I145, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I145), ta1_(t2), ta2_(v2) { }
};

class Task556 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task556(std::shared_ptr<TATensor<double,4>> I145, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I145), ta1_(t2), ta2_(v2) { }
};

class Task557 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task557(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> I148)
   : ta0_(I144), ta1_(Gamma49), ta2_(I148) { }
};

class Task558 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task558(std::shared_ptr<TATensor<double,4>> I148, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I148), ta1_(t2), ta2_(h1) { }
};

class Task559 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task559(std::shared_ptr<TATensor<double,4>> I148, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I148), ta1_(t2), ta2_(v2) { }
};

class Task560 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task560(std::shared_ptr<TATensor<double,4>> I148, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I148), ta1_(t2), ta2_(v2) { }
};

class Task561 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task561(std::shared_ptr<TATensor<double,4>> I148, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I148), ta1_(t2), ta2_(v2) { }
};

class Task562 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task562(std::shared_ptr<TATensor<double,4>> I148, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I148), ta1_(t2), ta2_(v2) { }
};

class Task563 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task563(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> I151)
   : ta0_(I144), ta1_(Gamma50), ta2_(I151) { }
};

class Task564 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task564(std::shared_ptr<TATensor<double,4>> I151, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I151), ta1_(t2), ta2_(h1) { }
};

class Task565 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task565(std::shared_ptr<TATensor<double,4>> I151, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I151), ta1_(t2), ta2_(h1) { }
};

class Task566 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task566(std::shared_ptr<TATensor<double,4>> I151, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1075)
   : ta0_(I151), ta1_(t2), ta2_(I1075) { }
};

class Task567 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task567(std::shared_ptr<TATensor<double,4>> I1075, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1075), ta1_(v2) { }
};

class Task568 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task568(std::shared_ptr<TATensor<double,4>> I151, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1078)
   : ta0_(I151), ta1_(t2), ta2_(I1078) { }
};

class Task569 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task569(std::shared_ptr<TATensor<double,4>> I1078, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1078), ta1_(v2) { }
};

class Task570 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task570(std::shared_ptr<TATensor<double,4>> I151, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I151), ta1_(t2), ta2_(v2) { }
};

class Task571 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task571(std::shared_ptr<TATensor<double,4>> I151, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I151), ta1_(t2), ta2_(v2) { }
};

class Task572 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task572(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,2>> I154)
   : ta0_(I144), ta1_(Gamma51), ta2_(I154) { }
};

class Task573 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task573(std::shared_ptr<TATensor<double,2>> I154, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I154), ta1_(t2), ta2_(h1) { }
};

class Task574 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task574(std::shared_ptr<TATensor<double,2>> I154, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I154), ta1_(t2), ta2_(h1) { }
};

class Task575 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task575(std::shared_ptr<TATensor<double,2>> I154, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I154), ta1_(t2), ta2_(v2) { }
};

class Task576 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task576(std::shared_ptr<TATensor<double,2>> I154, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I154), ta1_(t2), ta2_(v2) { }
};

class Task577 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task577(std::shared_ptr<TATensor<double,2>> I154, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I154), ta1_(t2), ta2_(v2) { }
};

class Task578 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task578(std::shared_ptr<TATensor<double,2>> I154, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I154), ta1_(t2), ta2_(v2) { }
};

class Task579 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task579(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I1026)
   : ta0_(I144), ta1_(v2), ta2_(I1026) { }
};

class Task580 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task580(std::shared_ptr<TATensor<double,6>> I1026, std::shared_ptr<TATensor<double,8>> Gamma339, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1026), ta1_(Gamma339), ta2_(t2) { }
};

class Task581 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task581(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma340, std::shared_ptr<TATensor<double,4>> I1029)
   : ta0_(I144), ta1_(Gamma340), ta2_(I1029) { }
};

class Task582 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task582(std::shared_ptr<TATensor<double,4>> I1029, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1029), ta1_(t2), ta2_(v2) { }
};

class Task583 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task583(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1032)
   : ta0_(I144), ta1_(t2), ta2_(I1032) { }
};

class Task584 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task584(std::shared_ptr<TATensor<double,6>> I1032, std::shared_ptr<TATensor<double,8>> Gamma341, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1032), ta1_(Gamma341), ta2_(v2) { }
};

class Task585 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task585(std::shared_ptr<TATensor<double,6>> I1032, std::shared_ptr<TATensor<double,8>> Gamma342, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1032), ta1_(Gamma342), ta2_(v2) { }
};

class Task586 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task586(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1044)
   : ta0_(I144), ta1_(t2), ta2_(I1044) { }
};

class Task587 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task587(std::shared_ptr<TATensor<double,6>> I1044, std::shared_ptr<TATensor<double,8>> Gamma345, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1044), ta1_(Gamma345), ta2_(v2) { }
};

class Task588 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task588(std::shared_ptr<TATensor<double,6>> I1044, std::shared_ptr<TATensor<double,8>> Gamma346, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1044), ta1_(Gamma346), ta2_(v2) { }
};

class Task589 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task589(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,8>> Gamma349, std::shared_ptr<TATensor<double,6>> I1056)
   : ta0_(I144), ta1_(Gamma349), ta2_(I1056) { }
};

class Task590 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task590(std::shared_ptr<TATensor<double,6>> I1056, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1057)
   : ta0_(I1056), ta1_(t2), ta2_(I1057) { }
};

class Task591 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task591(std::shared_ptr<TATensor<double,4>> I1057, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1057), ta1_(v2) { }
};

class Task592 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task592(std::shared_ptr<TATensor<double,6>> I1056, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1056), ta1_(t2), ta2_(v2) { }
};

class Task593 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task593(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,8>> Gamma350, std::shared_ptr<TATensor<double,6>> I1059)
   : ta0_(I144), ta1_(Gamma350), ta2_(I1059) { }
};

class Task594 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task594(std::shared_ptr<TATensor<double,6>> I1059, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1059), ta1_(t2), ta2_(v2) { }
};

class Task595 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task595(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,8>> Gamma351, std::shared_ptr<TATensor<double,6>> I1062)
   : ta0_(I144), ta1_(Gamma351), ta2_(I1062) { }
};

class Task596 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task596(std::shared_ptr<TATensor<double,6>> I1062, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1062), ta1_(t2), ta2_(v2) { }
};

class Task597 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task597(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma359, std::shared_ptr<TATensor<double,4>> I1086)
   : ta0_(I144), ta1_(Gamma359), ta2_(I1086) { }
};

class Task598 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task598(std::shared_ptr<TATensor<double,4>> I1086, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1086), ta1_(t2), ta2_(v2) { }
};

class Task599 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task599(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,8>> Gamma366, std::shared_ptr<TATensor<double,6>> I1107)
   : ta0_(I144), ta1_(Gamma366), ta2_(I1107) { }
};


}
}
}
#endif
#endif

