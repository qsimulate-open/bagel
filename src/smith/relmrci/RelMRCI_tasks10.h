//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks10.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS10_H
#define __SRC_SMITH_RelMRCI_TASKS10_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task450 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task450(std::shared_ptr<TATensor<std::complex<double>,2>> I127, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I127), ta1_(h1) { }
};

class Task451 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task451(std::shared_ptr<TATensor<std::complex<double>,2>> I127, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I877)
   : ta0_(I127), ta1_(Gamma27), ta2_(I877) { }
};

class Task452 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task452(std::shared_ptr<TATensor<std::complex<double>,4>> I877, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I877), ta1_(v2) { }
};

class Task453 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task453(std::shared_ptr<TATensor<std::complex<double>,2>> I127, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I127), ta1_(Gamma11), ta2_(v2) { }
};

class Task454 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task454(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I129)
   : ta0_(I108), ta1_(t2), ta2_(I129) { }
};

class Task455 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task455(std::shared_ptr<TATensor<std::complex<double>,2>> I129, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I129), ta1_(Gamma27), ta2_(h1) { }
};

class Task456 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task456(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I132)
   : ta0_(I108), ta1_(t2), ta2_(I132) { }
};

class Task457 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task457(std::shared_ptr<TATensor<std::complex<double>,2>> I132, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I132), ta1_(Gamma27), ta2_(h1) { }
};

class Task458 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task458(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I777)
   : ta0_(I108), ta1_(v2), ta2_(I777) { }
};

class Task459 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task459(std::shared_ptr<TATensor<std::complex<double>,2>> I777, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I777), ta1_(Gamma9), ta2_(t2) { }
};

class Task460 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task460(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I780)
   : ta0_(I108), ta1_(v2), ta2_(I780) { }
};

class Task461 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task461(std::shared_ptr<TATensor<std::complex<double>,2>> I780, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I780), ta1_(Gamma9), ta2_(t2) { }
};

class Task462 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task462(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I783)
   : ta0_(I108), ta1_(t2), ta2_(I783) { }
};

class Task463 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task463(std::shared_ptr<TATensor<std::complex<double>,2>> I783, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I783), ta1_(Gamma5), ta2_(v2) { }
};

class Task464 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task464(std::shared_ptr<TATensor<std::complex<double>,2>> I783, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I783), ta1_(Gamma160), ta2_(v2) { }
};

class Task465 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task465(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I786)
   : ta0_(I108), ta1_(t2), ta2_(I786) { }
};

class Task466 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task466(std::shared_ptr<TATensor<std::complex<double>,2>> I786, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I786), ta1_(Gamma5), ta2_(v2) { }
};

class Task467 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task467(std::shared_ptr<TATensor<std::complex<double>,2>> I786, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma160, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I786), ta1_(Gamma160), ta2_(v2) { }
};

class Task468 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task468(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I795)
   : ta0_(I108), ta1_(t2), ta2_(I795) { }
};

class Task469 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task469(std::shared_ptr<TATensor<std::complex<double>,4>> I795, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I796)
   : ta0_(I795), ta1_(Gamma11), ta2_(I796) { }
};

class Task470 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task470(std::shared_ptr<TATensor<std::complex<double>,4>> I796, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I796), ta1_(v2) { }
};

class Task471 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task471(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I798)
   : ta0_(I108), ta1_(t2), ta2_(I798) { }
};

class Task472 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task472(std::shared_ptr<TATensor<std::complex<double>,4>> I798, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I799)
   : ta0_(I798), ta1_(Gamma11), ta2_(I799) { }
};

class Task473 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task473(std::shared_ptr<TATensor<std::complex<double>,4>> I799, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I799), ta1_(v2) { }
};

class Task474 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task474(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I807)
   : ta0_(I108), ta1_(t2), ta2_(I807) { }
};

class Task475 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task475(std::shared_ptr<TATensor<std::complex<double>,4>> I807, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I808)
   : ta0_(I807), ta1_(Gamma11), ta2_(I808) { }
};

class Task476 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task476(std::shared_ptr<TATensor<std::complex<double>,4>> I808, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I808), ta1_(v2) { }
};

class Task477 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task477(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I810)
   : ta0_(I108), ta1_(t2), ta2_(I810) { }
};

class Task478 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task478(std::shared_ptr<TATensor<std::complex<double>,4>> I810, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I811)
   : ta0_(I810), ta1_(Gamma11), ta2_(I811) { }
};

class Task479 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task479(std::shared_ptr<TATensor<std::complex<double>,4>> I811, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I811), ta1_(v2) { }
};

class Task480 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task480(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I819)
   : ta0_(I108), ta1_(v2), ta2_(I819) { }
};

class Task481 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task481(std::shared_ptr<TATensor<std::complex<double>,4>> I819, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I819), ta1_(Gamma11), ta2_(t2) { }
};

class Task482 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task482(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I822)
   : ta0_(I108), ta1_(v2), ta2_(I822) { }
};

class Task483 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task483(std::shared_ptr<TATensor<std::complex<double>,4>> I822, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I822), ta1_(Gamma11), ta2_(t2) { }
};

class Task484 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task484(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I825)
   : ta0_(I108), ta1_(t2), ta2_(I825) { }
};

class Task485 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task485(std::shared_ptr<TATensor<std::complex<double>,4>> I825, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> I826)
   : ta0_(I825), ta1_(Gamma24), ta2_(I826) { }
};

class Task486 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task486(std::shared_ptr<TATensor<std::complex<double>,4>> I826, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I826), ta1_(v2) { }
};

class Task487 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task487(std::shared_ptr<TATensor<std::complex<double>,4>> I825, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I825), ta1_(Gamma9), ta2_(v2) { }
};

class Task488 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task488(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I828)
   : ta0_(I108), ta1_(t2), ta2_(I828) { }
};

class Task489 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task489(std::shared_ptr<TATensor<std::complex<double>,4>> I828, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> I829)
   : ta0_(I828), ta1_(Gamma24), ta2_(I829) { }
};

class Task490 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task490(std::shared_ptr<TATensor<std::complex<double>,4>> I829, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I829), ta1_(v2) { }
};

class Task491 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task491(std::shared_ptr<TATensor<std::complex<double>,4>> I828, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma9, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I828), ta1_(Gamma9), ta2_(v2) { }
};

class Task492 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task492(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I849)
   : ta0_(I108), ta1_(v2), ta2_(I849) { }
};

class Task493 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task493(std::shared_ptr<TATensor<std::complex<double>,2>> I849, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I849), ta1_(Gamma27), ta2_(t2) { }
};

class Task494 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task494(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I852)
   : ta0_(I108), ta1_(v2), ta2_(I852) { }
};

class Task495 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task495(std::shared_ptr<TATensor<std::complex<double>,2>> I852, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I852), ta1_(Gamma27), ta2_(t2) { }
};

class Task496 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task496(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I855)
   : ta0_(I108), ta1_(v2), ta2_(I855) { }
};

class Task497 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task497(std::shared_ptr<TATensor<std::complex<double>,2>> I855, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I855), ta1_(Gamma27), ta2_(t2) { }
};

class Task498 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task498(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I858)
   : ta0_(I108), ta1_(v2), ta2_(I858) { }
};

class Task499 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task499(std::shared_ptr<TATensor<std::complex<double>,2>> I858, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I858), ta1_(Gamma27), ta2_(t2) { }
};


}
}
}
#endif
#endif

