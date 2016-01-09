//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks7.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS7_H
#define __SRC_SMITH_RelMRCI_TASKS7_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task300 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task300(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I67), ta1_(t2), ta2_(h1) { }
};

class Task301 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task301(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I67), ta1_(t2), ta2_(h1) { }
};

class Task302 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task302(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task303 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task303(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task304 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task304(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task305 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task305(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task306 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task306(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task307 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task307(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task308 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task308(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task309 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task309(std::shared_ptr<TATensor<std::complex<double>,4>> I67, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I67), ta1_(t2), ta2_(v2) { }
};

class Task310 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task310(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> I73)
   : ta0_(I63), ta1_(Gamma24), ta2_(I73) { }
};

class Task311 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task311(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I73), ta1_(t2), ta2_(h1) { }
};

class Task312 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task312(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I73), ta1_(t2), ta2_(h1) { }
};

class Task313 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task313(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I73), ta1_(t2), ta2_(h1) { }
};

class Task314 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task314(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I73), ta1_(t2), ta2_(h1) { }
};

class Task315 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task315(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I73), ta1_(t2), ta2_(v2) { }
};

class Task316 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task316(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I73), ta1_(t2), ta2_(v2) { }
};

class Task317 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task317(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I631)
   : ta0_(I73), ta1_(t2), ta2_(I631) { }
};

class Task318 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task318(std::shared_ptr<TATensor<std::complex<double>,4>> I631, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I631), ta1_(v2) { }
};

class Task319 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task319(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I634)
   : ta0_(I73), ta1_(t2), ta2_(I634) { }
};

class Task320 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task320(std::shared_ptr<TATensor<std::complex<double>,4>> I634, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I634), ta1_(v2) { }
};

class Task321 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task321(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I73), ta1_(t2), ta2_(v2) { }
};

class Task322 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task322(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I73), ta1_(t2), ta2_(v2) { }
};

class Task323 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task323(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I679)
   : ta0_(I73), ta1_(t2), ta2_(I679) { }
};

class Task324 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task324(std::shared_ptr<TATensor<std::complex<double>,4>> I679, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I679), ta1_(v2) { }
};

class Task325 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task325(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I682)
   : ta0_(I73), ta1_(t2), ta2_(I682) { }
};

class Task326 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task326(std::shared_ptr<TATensor<std::complex<double>,4>> I682, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I682), ta1_(v2) { }
};

class Task327 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task327(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I73), ta1_(t2), ta2_(v2) { }
};

class Task328 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task328(std::shared_ptr<TATensor<std::complex<double>,4>> I73, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I73), ta1_(t2), ta2_(v2) { }
};

class Task329 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task329(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I79)
   : ta0_(I63), ta1_(h1), ta2_(I79) { }
};

class Task330 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task330(std::shared_ptr<TATensor<std::complex<double>,4>> I79, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma26, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I79), ta1_(Gamma26), ta2_(t2) { }
};

class Task331 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task331(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> I82)
   : ta0_(I63), ta1_(Gamma27), ta2_(I82) { }
};

class Task332 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task332(std::shared_ptr<TATensor<std::complex<double>,2>> I82, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I82), ta1_(t2), ta2_(h1) { }
};

class Task333 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task333(std::shared_ptr<TATensor<std::complex<double>,2>> I82, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I82), ta1_(t2), ta2_(h1) { }
};

class Task334 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task334(std::shared_ptr<TATensor<std::complex<double>,2>> I82, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task335 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task335(std::shared_ptr<TATensor<std::complex<double>,2>> I82, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task336 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task336(std::shared_ptr<TATensor<std::complex<double>,2>> I82, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task337 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task337(std::shared_ptr<TATensor<std::complex<double>,2>> I82, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task338 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task338(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,6>> I543)
   : ta0_(I63), ta1_(v2), ta2_(I543) { }
};

class Task339 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task339(std::shared_ptr<TATensor<std::complex<double>,6>> I543, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma84, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I543), ta1_(Gamma84), ta2_(t2) { }
};

class Task340 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task340(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,6>> I546)
   : ta0_(I63), ta1_(v2), ta2_(I546) { }
};

class Task341 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task341(std::shared_ptr<TATensor<std::complex<double>,6>> I546, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma179, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I546), ta1_(Gamma179), ta2_(t2) { }
};

class Task342 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task342(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,6>> I549)
   : ta0_(I63), ta1_(v2), ta2_(I549) { }
};

class Task343 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task343(std::shared_ptr<TATensor<std::complex<double>,6>> I549, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma77, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I549), ta1_(Gamma77), ta2_(t2) { }
};

class Task344 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task344(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I552)
   : ta0_(I63), ta1_(v2), ta2_(I552) { }
};

class Task345 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task345(std::shared_ptr<TATensor<std::complex<double>,4>> I552, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I552), ta1_(Gamma4), ta2_(t2) { }
};

class Task346 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task346(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I555)
   : ta0_(I63), ta1_(v2), ta2_(I555) { }
};

class Task347 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task347(std::shared_ptr<TATensor<std::complex<double>,4>> I555, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I555), ta1_(Gamma4), ta2_(t2) { }
};

class Task348 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task348(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I558)
   : ta0_(I63), ta1_(t2), ta2_(I558) { }
};

class Task349 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task349(std::shared_ptr<TATensor<std::complex<double>,4>> I558, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma183, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I558), ta1_(Gamma183), ta2_(v2) { }
};


}
}
}
#endif
#endif

