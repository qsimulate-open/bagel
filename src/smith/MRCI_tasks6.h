//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks6.h
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

#ifndef __SRC_SMITH_MRCI_TASKS6_H
#define __SRC_SMITH_MRCI_TASKS6_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task250 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task250(std::shared_ptr<TATensor<double,4>> I399, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I399), ta1_(Gamma0), ta2_(t2) { }
};

class Task251 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task251(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I402)
   : ta0_(I27), ta1_(v2), ta2_(I402) { }
};

class Task252 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task252(std::shared_ptr<TATensor<double,4>> I402, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I402), ta1_(Gamma2), ta2_(t2) { }
};

class Task253 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task253(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I405)
   : ta0_(I27), ta1_(v2), ta2_(I405) { }
};

class Task254 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task254(std::shared_ptr<TATensor<double,4>> I405, std::shared_ptr<TATensor<double,6>> Gamma132, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I405), ta1_(Gamma132), ta2_(t2) { }
};

class Task255 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task255(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I408)
   : ta0_(I27), ta1_(v2), ta2_(I408) { }
};

class Task256 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task256(std::shared_ptr<TATensor<double,4>> I408, std::shared_ptr<TATensor<double,6>> Gamma132, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I408), ta1_(Gamma132), ta2_(t2) { }
};

class Task257 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task257(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I411)
   : ta0_(I27), ta1_(v2), ta2_(I411) { }
};

class Task258 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task258(std::shared_ptr<TATensor<double,4>> I411, std::shared_ptr<TATensor<double,6>> Gamma1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I411), ta1_(Gamma1), ta2_(t2) { }
};

class Task259 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task259(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I414)
   : ta0_(I27), ta1_(v2), ta2_(I414) { }
};

class Task260 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task260(std::shared_ptr<TATensor<double,4>> I414, std::shared_ptr<TATensor<double,6>> Gamma87, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I414), ta1_(Gamma87), ta2_(t2) { }
};

class Task261 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task261(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I417)
   : ta0_(I27), ta1_(v2), ta2_(I417) { }
};

class Task262 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task262(std::shared_ptr<TATensor<double,4>> I417, std::shared_ptr<TATensor<double,6>> Gamma132, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I417), ta1_(Gamma132), ta2_(t2) { }
};

class Task263 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task263(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I420)
   : ta0_(I27), ta1_(v2), ta2_(I420) { }
};

class Task264 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task264(std::shared_ptr<TATensor<double,4>> I420, std::shared_ptr<TATensor<double,6>> Gamma137, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I420), ta1_(Gamma137), ta2_(t2) { }
};

class Task265 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task265(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I423)
   : ta0_(I27), ta1_(v2), ta2_(I423) { }
};

class Task266 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task266(std::shared_ptr<TATensor<double,4>> I423, std::shared_ptr<TATensor<double,6>> Gamma132, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I423), ta1_(Gamma132), ta2_(t2) { }
};

class Task267 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task267(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I426)
   : ta0_(I27), ta1_(v2), ta2_(I426) { }
};

class Task268 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task268(std::shared_ptr<TATensor<double,4>> I426, std::shared_ptr<TATensor<double,6>> Gamma132, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I426), ta1_(Gamma132), ta2_(t2) { }
};

class Task269 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task269(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I429)
   : ta0_(I27), ta1_(v2), ta2_(I429) { }
};

class Task270 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task270(std::shared_ptr<TATensor<double,2>> I429, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I429), ta1_(Gamma10), ta2_(t2) { }
};

class Task271 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task271(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I432)
   : ta0_(I27), ta1_(v2), ta2_(I432) { }
};

class Task272 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task272(std::shared_ptr<TATensor<double,2>> I432, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I432), ta1_(Gamma10), ta2_(t2) { }
};

class Task273 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task273(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I435)
   : ta0_(I27), ta1_(t2), ta2_(I435) { }
};

class Task274 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task274(std::shared_ptr<TATensor<double,4>> I435, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> I436)
   : ta0_(I435), ta1_(Gamma197), ta2_(I436) { }
};

class Task275 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task275(std::shared_ptr<TATensor<double,4>> I436, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I436), ta1_(v2) { }
};

class Task276 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task276(std::shared_ptr<TATensor<double,4>> I435, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I435), ta1_(Gamma0), ta2_(v2) { }
};

class Task277 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task277(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I438)
   : ta0_(I27), ta1_(t2), ta2_(I438) { }
};

class Task278 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task278(std::shared_ptr<TATensor<double,4>> I438, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> I439)
   : ta0_(I438), ta1_(Gamma197), ta2_(I439) { }
};

class Task279 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task279(std::shared_ptr<TATensor<double,4>> I439, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I439), ta1_(v2) { }
};

class Task280 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task280(std::shared_ptr<TATensor<double,4>> I438, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I438), ta1_(Gamma2), ta2_(v2) { }
};

class Task281 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task281(std::shared_ptr<TATensor<double,4>> I438, std::shared_ptr<TATensor<double,4>> Gamma155, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I438), ta1_(Gamma155), ta2_(v2) { }
};

class Task282 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task282(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I441)
   : ta0_(I27), ta1_(t2), ta2_(I441) { }
};

class Task283 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task283(std::shared_ptr<TATensor<double,4>> I441, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> I442)
   : ta0_(I441), ta1_(Gamma197), ta2_(I442) { }
};

class Task284 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task284(std::shared_ptr<TATensor<double,4>> I442, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I442), ta1_(v2) { }
};

class Task285 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task285(std::shared_ptr<TATensor<double,4>> I441, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I441), ta1_(Gamma2), ta2_(v2) { }
};

class Task286 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task286(std::shared_ptr<TATensor<double,4>> I441, std::shared_ptr<TATensor<double,4>> Gamma155, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I441), ta1_(Gamma155), ta2_(v2) { }
};

class Task287 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task287(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I444)
   : ta0_(I27), ta1_(t2), ta2_(I444) { }
};

class Task288 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task288(std::shared_ptr<TATensor<double,4>> I444, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> I445)
   : ta0_(I444), ta1_(Gamma197), ta2_(I445) { }
};

class Task289 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task289(std::shared_ptr<TATensor<double,4>> I445, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I445), ta1_(v2) { }
};

class Task290 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task290(std::shared_ptr<TATensor<double,4>> I444, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I444), ta1_(Gamma2), ta2_(v2) { }
};

class Task291 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task291(std::shared_ptr<TATensor<double,4>> I444, std::shared_ptr<TATensor<double,4>> Gamma155, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I444), ta1_(Gamma155), ta2_(v2) { }
};

class Task292 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task292(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I447)
   : ta0_(I27), ta1_(t2), ta2_(I447) { }
};

class Task293 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task293(std::shared_ptr<TATensor<double,4>> I447, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> I448)
   : ta0_(I447), ta1_(Gamma197), ta2_(I448) { }
};

class Task294 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task294(std::shared_ptr<TATensor<double,4>> I448, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I448), ta1_(v2) { }
};

class Task295 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task295(std::shared_ptr<TATensor<double,4>> I447, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I447), ta1_(Gamma2), ta2_(v2) { }
};

class Task296 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task296(std::shared_ptr<TATensor<double,4>> I447, std::shared_ptr<TATensor<double,4>> Gamma155, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I447), ta1_(Gamma155), ta2_(v2) { }
};

class Task297 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task297(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I450)
   : ta0_(I27), ta1_(t2), ta2_(I450) { }
};

class Task298 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task298(std::shared_ptr<TATensor<double,4>> I450, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> I451)
   : ta0_(I450), ta1_(Gamma197), ta2_(I451) { }
};

class Task299 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task299(std::shared_ptr<TATensor<double,4>> I451, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I451), ta1_(v2) { }
};


}
}
}
#endif
#endif

