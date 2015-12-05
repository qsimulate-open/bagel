//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks12.h
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
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task550(std::shared_ptr<TATensor<double,4>> I763, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I763), ta1_(t2) { }
};

class Task551 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task551(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I766)
   : ta0_(proj), ta1_(I766) { }
};

class Task552 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task552(std::shared_ptr<TATensor<double,4>> I766, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I766), ta1_(Gamma60), ta2_(t2) { }
};

class Task553 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task553(std::shared_ptr<TATensor<double,1>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task554 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,1>> ta1_;
    void compute_() override;
  public:
    Task554(std::shared_ptr<TATensor<double,1>> proj, std::shared_ptr<TATensor<double,1>> I768)
   : ta0_(proj), ta1_(I768) { }
};

class Task555 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task555(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma270, std::shared_ptr<TATensor<double,4>> I769)
   : ta0_(I768), ta1_(Gamma270), ta2_(I769) { }
};

class Task556 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task556(std::shared_ptr<TATensor<double,4>> I769, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I769), ta1_(t2) { }
};

class Task557 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task557(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma271, std::shared_ptr<TATensor<double,4>> I772)
   : ta0_(I768), ta1_(Gamma271), ta2_(I772) { }
};

class Task558 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task558(std::shared_ptr<TATensor<double,4>> I772, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I773)
   : ta0_(I772), ta1_(t2), ta2_(I773) { }
};

class Task559 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task559(std::shared_ptr<TATensor<double,4>> I773, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I773), ta1_(f1), ta2_(t2) { }
};

class Task560 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task560(std::shared_ptr<TATensor<double,4>> I772, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I772), ta1_(t2), e0_(e) { }
};

class Task561 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task561(std::shared_ptr<TATensor<double,4>> I772, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I772), ta1_(v2), ta2_(t2) { }
};

class Task562 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task562(std::shared_ptr<TATensor<double,4>> I772, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I772), ta1_(v2), ta2_(t2) { }
};

class Task563 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task563(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma272, std::shared_ptr<TATensor<double,6>> I776)
   : ta0_(I768), ta1_(Gamma272), ta2_(I776) { }
};

class Task564 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task564(std::shared_ptr<TATensor<double,6>> I776, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I777)
   : ta0_(I776), ta1_(t2), ta2_(I777) { }
};

class Task565 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task565(std::shared_ptr<TATensor<double,4>> I777, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I777), ta1_(f1), ta2_(t2) { }
};

class Task566 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task566(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma273, std::shared_ptr<TATensor<double,4>> I780)
   : ta0_(I768), ta1_(Gamma273), ta2_(I780) { }
};

class Task567 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task567(std::shared_ptr<TATensor<double,4>> I780, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I781)
   : ta0_(I780), ta1_(t2), ta2_(I781) { }
};

class Task568 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task568(std::shared_ptr<TATensor<double,4>> I781, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I781), ta1_(t2), ta2_(f1) { }
};

class Task569 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task569(std::shared_ptr<TATensor<double,4>> I780, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I812)
   : ta0_(I780), ta1_(t2), ta2_(I812) { }
};

class Task570 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task570(std::shared_ptr<TATensor<double,4>> I812, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I812), ta1_(f1), ta2_(t2) { }
};

class Task571 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task571(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma274, std::shared_ptr<TATensor<double,6>> I784)
   : ta0_(I768), ta1_(Gamma274), ta2_(I784) { }
};

class Task572 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task572(std::shared_ptr<TATensor<double,6>> I784, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I785)
   : ta0_(I784), ta1_(t2), ta2_(I785) { }
};

class Task573 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task573(std::shared_ptr<TATensor<double,4>> I785, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I785), ta1_(t2), ta2_(f1) { }
};

class Task574 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task574(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma275, std::shared_ptr<TATensor<double,6>> I788)
   : ta0_(I768), ta1_(Gamma275), ta2_(I788) { }
};

class Task575 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task575(std::shared_ptr<TATensor<double,6>> I788, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I788), ta1_(t2) { }
};

class Task576 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task576(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma276, std::shared_ptr<TATensor<double,6>> I791)
   : ta0_(I768), ta1_(Gamma276), ta2_(I791) { }
};

class Task577 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task577(std::shared_ptr<TATensor<double,6>> I791, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I792)
   : ta0_(I791), ta1_(t2), ta2_(I792) { }
};

class Task578 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task578(std::shared_ptr<TATensor<double,4>> I792, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I792), ta1_(f1), ta2_(t2) { }
};

class Task579 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task579(std::shared_ptr<TATensor<double,6>> I791, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I808)
   : ta0_(I791), ta1_(t2), ta2_(I808) { }
};

class Task580 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task580(std::shared_ptr<TATensor<double,4>> I808, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I808), ta1_(t2), ta2_(f1) { }
};

class Task581 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task581(std::shared_ptr<TATensor<double,6>> I791, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I932)
   : ta0_(I791), ta1_(t2), ta2_(I932) { }
};

class Task582 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task582(std::shared_ptr<TATensor<double,4>> I932, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I932), ta1_(t2), e0_(e) { }
};

class Task583 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task583(std::shared_ptr<TATensor<double,4>> I932, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I932), ta1_(f1), ta2_(t2) { }
};

class Task584 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task584(std::shared_ptr<TATensor<double,6>> I791, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I791), ta1_(v2), ta2_(t2) { }
};

class Task585 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task585(std::shared_ptr<TATensor<double,6>> I791, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I791), ta1_(v2), ta2_(t2) { }
};

class Task586 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task586(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma277, std::shared_ptr<TATensor<double,4>> I795)
   : ta0_(I768), ta1_(Gamma277), ta2_(I795) { }
};

class Task587 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task587(std::shared_ptr<TATensor<double,4>> I795, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I796)
   : ta0_(I795), ta1_(t2), ta2_(I796) { }
};

class Task588 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task588(std::shared_ptr<TATensor<double,2>> I796, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I796), ta1_(t2), ta2_(f1) { }
};

class Task589 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task589(std::shared_ptr<TATensor<double,2>> I796, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I796), ta1_(t2), ta2_(f1) { }
};

class Task590 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task590(std::shared_ptr<TATensor<double,4>> I795, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I886)
   : ta0_(I795), ta1_(t2), ta2_(I886) { }
};

class Task591 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task591(std::shared_ptr<TATensor<double,4>> I886, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I886), ta1_(t2), ta2_(f1) { }
};

class Task592 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task592(std::shared_ptr<TATensor<double,4>> I795, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I936)
   : ta0_(I795), ta1_(t2), ta2_(I936) { }
};

class Task593 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task593(std::shared_ptr<TATensor<double,4>> I936, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I936), ta1_(t2), ta2_(f1) { }
};

class Task594 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task594(std::shared_ptr<TATensor<double,4>> I936, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I936), ta1_(t2), ta2_(f1) { }
};

class Task595 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task595(std::shared_ptr<TATensor<double,4>> I795, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I795), ta1_(v2), ta2_(t2) { }
};

class Task596 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task596(std::shared_ptr<TATensor<double,4>> I795, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I795), ta1_(h1), ta2_(t2) { }
};

class Task597 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task597(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma279, std::shared_ptr<TATensor<double,6>> I803)
   : ta0_(I768), ta1_(Gamma279), ta2_(I803) { }
};

class Task598 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task598(std::shared_ptr<TATensor<double,6>> I803, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I804)
   : ta0_(I803), ta1_(t2), ta2_(I804) { }
};

class Task599 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task599(std::shared_ptr<TATensor<double,4>> I804, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I804), ta1_(t2), ta2_(f1) { }
};


}
}
}
#endif
#endif

