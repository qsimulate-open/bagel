//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks9.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS9_H
#define __SRC_SMITH_RelMRCI_TASKS9_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task400 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task400(std::shared_ptr<TATensor<std::complex<double>,2>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I100), ta1_(t2), ta2_(v2) { }
};

class Task401 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task401(std::shared_ptr<TATensor<std::complex<double>,2>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I100), ta1_(t2), ta2_(v2) { }
};

class Task402 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task402(std::shared_ptr<TATensor<std::complex<double>,2>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I100), ta1_(t2), ta2_(v2) { }
};

class Task403 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task403(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,6>> I699)
   : ta0_(I93), ta1_(v2), ta2_(I699) { }
};

class Task404 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task404(std::shared_ptr<TATensor<std::complex<double>,6>> I699, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma230, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I699), ta1_(Gamma230), ta2_(t2) { }
};

class Task405 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task405(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma231, std::shared_ptr<TATensor<std::complex<double>,4>> I702)
   : ta0_(I93), ta1_(Gamma231), ta2_(I702) { }
};

class Task406 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task406(std::shared_ptr<TATensor<std::complex<double>,4>> I702, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I702), ta1_(t2), ta2_(v2) { }
};

class Task407 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task407(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I705)
   : ta0_(I93), ta1_(t2), ta2_(I705) { }
};

class Task408 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task408(std::shared_ptr<TATensor<std::complex<double>,6>> I705, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma232, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I705), ta1_(Gamma232), ta2_(v2) { }
};

class Task409 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task409(std::shared_ptr<TATensor<std::complex<double>,6>> I705, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma233, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I705), ta1_(Gamma233), ta2_(v2) { }
};

class Task410 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task410(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma236, std::shared_ptr<TATensor<std::complex<double>,6>> I717)
   : ta0_(I93), ta1_(Gamma236), ta2_(I717) { }
};

class Task411 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task411(std::shared_ptr<TATensor<std::complex<double>,6>> I717, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I718)
   : ta0_(I717), ta1_(t2), ta2_(I718) { }
};

class Task412 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task412(std::shared_ptr<TATensor<std::complex<double>,4>> I718, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I718), ta1_(v2) { }
};

class Task413 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task413(std::shared_ptr<TATensor<std::complex<double>,6>> I717, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I717), ta1_(t2), ta2_(v2) { }
};

class Task414 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task414(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma237, std::shared_ptr<TATensor<std::complex<double>,6>> I720)
   : ta0_(I93), ta1_(Gamma237), ta2_(I720) { }
};

class Task415 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task415(std::shared_ptr<TATensor<std::complex<double>,6>> I720, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I720), ta1_(t2), ta2_(v2) { }
};

class Task416 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task416(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma238, std::shared_ptr<TATensor<std::complex<double>,6>> I723)
   : ta0_(I93), ta1_(Gamma238), ta2_(I723) { }
};

class Task417 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task417(std::shared_ptr<TATensor<std::complex<double>,6>> I723, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I723), ta1_(t2), ta2_(v2) { }
};

class Task418 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task418(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma245, std::shared_ptr<TATensor<std::complex<double>,4>> I744)
   : ta0_(I93), ta1_(Gamma245), ta2_(I744) { }
};

class Task419 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task419(std::shared_ptr<TATensor<std::complex<double>,4>> I744, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I744), ta1_(t2), ta2_(v2) { }
};

class Task420 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task420(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma246, std::shared_ptr<TATensor<std::complex<double>,4>> I747)
   : ta0_(I93), ta1_(Gamma246), ta2_(I747) { }
};

class Task421 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task421(std::shared_ptr<TATensor<std::complex<double>,4>> I747, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I747), ta1_(t2), ta2_(v2) { }
};

class Task422 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task422(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma253, std::shared_ptr<TATensor<std::complex<double>,6>> I768)
   : ta0_(I93), ta1_(Gamma253), ta2_(I768) { }
};

class Task423 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task423(std::shared_ptr<TATensor<std::complex<double>,6>> I768, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I768), ta1_(t2), ta2_(v2) { }
};

class Task424 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task424(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma422, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I93), ta1_(Gamma422), ta2_(t2) { }
};

class Task425 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task425(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma423, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I93), ta1_(Gamma423), ta2_(t2) { }
};

class Task426 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task426(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I108)
   : ta0_(proj), ta1_(I108) { }
};

class Task427 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task427(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I109)
   : ta0_(I108), ta1_(t2), ta2_(I109) { }
};

class Task428 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task428(std::shared_ptr<TATensor<std::complex<double>,2>> I109, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I109), ta1_(Gamma11), ta2_(h1) { }
};

class Task429 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task429(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I112)
   : ta0_(I108), ta1_(t2), ta2_(I112) { }
};

class Task430 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task430(std::shared_ptr<TATensor<std::complex<double>,2>> I112, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I112), ta1_(Gamma11), ta2_(h1) { }
};

class Task431 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task431(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,2>> I115)
   : ta0_(I108), ta1_(h1), ta2_(I115) { }
};

class Task432 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task432(std::shared_ptr<TATensor<std::complex<double>,2>> I115, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I115), ta1_(Gamma27), ta2_(t2) { }
};

class Task433 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task433(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,2>> I118)
   : ta0_(I108), ta1_(h1), ta2_(I118) { }
};

class Task434 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task434(std::shared_ptr<TATensor<std::complex<double>,2>> I118, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I118), ta1_(Gamma27), ta2_(t2) { }
};

class Task435 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task435(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I121)
   : ta0_(I108), ta1_(t2), ta2_(I121) { }
};

class Task436 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task436(std::shared_ptr<TATensor<std::complex<double>,2>> I121, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I121), ta1_(h1) { }
};

class Task437 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task437(std::shared_ptr<TATensor<std::complex<double>,2>> I121, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I868)
   : ta0_(I121), ta1_(Gamma27), ta2_(I868) { }
};

class Task438 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task438(std::shared_ptr<TATensor<std::complex<double>,4>> I868, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I868), ta1_(v2) { }
};

class Task439 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task439(std::shared_ptr<TATensor<std::complex<double>,2>> I121, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I121), ta1_(Gamma11), ta2_(v2) { }
};

class Task440 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task440(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I123)
   : ta0_(I108), ta1_(t2), ta2_(I123) { }
};

class Task441 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task441(std::shared_ptr<TATensor<std::complex<double>,2>> I123, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I123), ta1_(h1) { }
};

class Task442 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task442(std::shared_ptr<TATensor<std::complex<double>,2>> I123, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I871)
   : ta0_(I123), ta1_(Gamma27), ta2_(I871) { }
};

class Task443 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task443(std::shared_ptr<TATensor<std::complex<double>,4>> I871, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I871), ta1_(v2) { }
};

class Task444 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task444(std::shared_ptr<TATensor<std::complex<double>,2>> I123, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I123), ta1_(Gamma11), ta2_(v2) { }
};

class Task445 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task445(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I125)
   : ta0_(I108), ta1_(t2), ta2_(I125) { }
};

class Task446 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    void compute_() override;
  public:
    Task446(std::shared_ptr<TATensor<std::complex<double>,2>> I125, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I125), ta1_(h1) { }
};

class Task447 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task447(std::shared_ptr<TATensor<std::complex<double>,2>> I125, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I874)
   : ta0_(I125), ta1_(Gamma27), ta2_(I874) { }
};

class Task448 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task448(std::shared_ptr<TATensor<std::complex<double>,4>> I874, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I874), ta1_(v2) { }
};

class Task449 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task449(std::shared_ptr<TATensor<std::complex<double>,2>> I125, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I125), ta1_(Gamma11), ta2_(v2) { }
};


}
}
}
#endif
#endif

