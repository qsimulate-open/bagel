//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks10.h
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

#ifndef __SRC_SMITH_MRCI_TASKS10_H
#define __SRC_SMITH_MRCI_TASKS10_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task450 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task450(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma572, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I72), ta1_(Gamma572), ta2_(t2) { }
};

class Task451 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task451(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma573, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I72), ta1_(Gamma573), ta2_(t2) { }
};

class Task452 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task452(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I108)
   : ta0_(proj), ta1_(I108) { }
};

class Task453 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task453(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I109)
   : ta0_(I108), ta1_(h1), ta2_(I109) { }
};

class Task454 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task454(std::shared_ptr<TATensor<double,4>> I109, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I109), ta1_(Gamma4), ta2_(t2) { }
};

class Task455 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task455(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,4>> I112)
   : ta0_(I108), ta1_(Gamma5), ta2_(I112) { }
};

class Task456 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task456(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I112), ta1_(t2), ta2_(h1) { }
};

class Task457 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task457(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I112), ta1_(t2), ta2_(h1) { }
};

class Task458 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task458(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task459 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task459(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task460 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task460(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task461 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task461(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task462 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task462(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task463 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task463(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task464 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task464(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task465 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task465(std::shared_ptr<TATensor<double,4>> I112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I112), ta1_(t2), ta2_(v2) { }
};

class Task466 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task466(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I118)
   : ta0_(I108), ta1_(Gamma29), ta2_(I118) { }
};

class Task467 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task467(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I118), ta1_(t2), ta2_(h1) { }
};

class Task468 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task468(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I118), ta1_(t2), ta2_(h1) { }
};

class Task469 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task469(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I118), ta1_(t2), ta2_(h1) { }
};

class Task470 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task470(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I118), ta1_(t2), ta2_(h1) { }
};

class Task471 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task471(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I118), ta1_(t2), ta2_(h1) { }
};

class Task472 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task472(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I118), ta1_(t2), ta2_(h1) { }
};

class Task473 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task473(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task474 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task474(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task475 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task475(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task476 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task476(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task477 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task477(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I958)
   : ta0_(I118), ta1_(t2), ta2_(I958) { }
};

class Task478 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task478(std::shared_ptr<TATensor<double,4>> I958, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I958), ta1_(v2) { }
};

class Task479 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task479(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I961)
   : ta0_(I118), ta1_(t2), ta2_(I961) { }
};

class Task480 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task480(std::shared_ptr<TATensor<double,4>> I961, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I961), ta1_(v2) { }
};

class Task481 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task481(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task482 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task482(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task483 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task483(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1006)
   : ta0_(I118), ta1_(t2), ta2_(I1006) { }
};

class Task484 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task484(std::shared_ptr<TATensor<double,4>> I1006, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1006), ta1_(v2) { }
};

class Task485 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task485(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1009)
   : ta0_(I118), ta1_(t2), ta2_(I1009) { }
};

class Task486 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task486(std::shared_ptr<TATensor<double,4>> I1009, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1009), ta1_(v2) { }
};

class Task487 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task487(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task488 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task488(std::shared_ptr<TATensor<double,4>> I118, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I118), ta1_(t2), ta2_(v2) { }
};

class Task489 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task489(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I130)
   : ta0_(I108), ta1_(h1), ta2_(I130) { }
};

class Task490 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task490(std::shared_ptr<TATensor<double,4>> I130, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I130), ta1_(Gamma252), ta2_(t2) { }
};

class Task491 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task491(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> I133)
   : ta0_(I108), ta1_(Gamma32), ta2_(I133) { }
};

class Task492 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task492(std::shared_ptr<TATensor<double,2>> I133, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I133), ta1_(t2), ta2_(h1) { }
};

class Task493 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task493(std::shared_ptr<TATensor<double,2>> I133, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I133), ta1_(t2), ta2_(h1) { }
};

class Task494 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task494(std::shared_ptr<TATensor<double,2>> I133, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I133), ta1_(t2), ta2_(v2) { }
};

class Task495 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task495(std::shared_ptr<TATensor<double,2>> I133, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I133), ta1_(t2), ta2_(v2) { }
};

class Task496 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task496(std::shared_ptr<TATensor<double,2>> I133, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I133), ta1_(t2), ta2_(v2) { }
};

class Task497 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task497(std::shared_ptr<TATensor<double,2>> I133, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I133), ta1_(t2), ta2_(v2) { }
};

class Task498 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task498(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I840)
   : ta0_(I108), ta1_(v2), ta2_(I840) { }
};

class Task499 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task499(std::shared_ptr<TATensor<double,6>> I840, std::shared_ptr<TATensor<double,6>> Gamma107, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I840), ta1_(Gamma107), ta2_(t2) { }
};


}
}
}
#endif
#endif

