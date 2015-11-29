//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks9.h
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

#ifndef __SRC_SMITH_MRCI_TASKS9_H
#define __SRC_SMITH_MRCI_TASKS9_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task400 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task400(std::shared_ptr<TATensor<double,4>> I663, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I663), ta1_(Gamma4), ta2_(t2) { }
};

class Task401 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task401(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I666)
   : ta0_(I72), ta1_(v2), ta2_(I666) { }
};

class Task402 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task402(std::shared_ptr<TATensor<double,4>> I666, std::shared_ptr<TATensor<double,6>> Gamma24, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I666), ta1_(Gamma24), ta2_(t2) { }
};

class Task403 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task403(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I669)
   : ta0_(I72), ta1_(t2), ta2_(I669) { }
};

class Task404 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task404(std::shared_ptr<TATensor<double,4>> I669, std::shared_ptr<TATensor<double,6>> Gamma220, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I669), ta1_(Gamma220), ta2_(v2) { }
};

class Task405 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task405(std::shared_ptr<TATensor<double,4>> I669, std::shared_ptr<TATensor<double,6>> Gamma222, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I669), ta1_(Gamma222), ta2_(v2) { }
};

class Task406 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task406(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I672)
   : ta0_(I72), ta1_(t2), ta2_(I672) { }
};

class Task407 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task407(std::shared_ptr<TATensor<double,4>> I672, std::shared_ptr<TATensor<double,6>> Gamma221, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I672), ta1_(Gamma221), ta2_(v2) { }
};

class Task408 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task408(std::shared_ptr<TATensor<double,4>> I672, std::shared_ptr<TATensor<double,6>> Gamma104, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I672), ta1_(Gamma104), ta2_(v2) { }
};

class Task409 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task409(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I699)
   : ta0_(I72), ta1_(t2), ta2_(I699) { }
};

class Task410 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task410(std::shared_ptr<TATensor<double,6>> I699, std::shared_ptr<TATensor<double,6>> Gamma230, std::shared_ptr<TATensor<double,4>> I700)
   : ta0_(I699), ta1_(Gamma230), ta2_(I700) { }
};

class Task411 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task411(std::shared_ptr<TATensor<double,4>> I700, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I700), ta1_(v2) { }
};

class Task412 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task412(std::shared_ptr<TATensor<double,6>> I699, std::shared_ptr<TATensor<double,6>> Gamma232, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I699), ta1_(Gamma232), ta2_(v2) { }
};

class Task413 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task413(std::shared_ptr<TATensor<double,6>> I699, std::shared_ptr<TATensor<double,6>> Gamma234, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I699), ta1_(Gamma234), ta2_(v2) { }
};

class Task414 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task414(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma230, std::shared_ptr<TATensor<double,6>> I702)
   : ta0_(I72), ta1_(Gamma230), ta2_(I702) { }
};

class Task415 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task415(std::shared_ptr<TATensor<double,6>> I702, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I703)
   : ta0_(I702), ta1_(t2), ta2_(I703) { }
};

class Task416 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task416(std::shared_ptr<TATensor<double,4>> I703, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I703), ta1_(v2) { }
};

class Task417 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task417(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma233, std::shared_ptr<TATensor<double,6>> I708)
   : ta0_(I72), ta1_(Gamma233), ta2_(I708) { }
};

class Task418 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task418(std::shared_ptr<TATensor<double,6>> I708, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I708), ta1_(t2), ta2_(v2) { }
};

class Task419 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task419(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma235, std::shared_ptr<TATensor<double,6>> I714)
   : ta0_(I72), ta1_(Gamma235), ta2_(I714) { }
};

class Task420 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task420(std::shared_ptr<TATensor<double,6>> I714, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I714), ta1_(t2), ta2_(v2) { }
};

class Task421 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task421(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I729)
   : ta0_(I72), ta1_(t2), ta2_(I729) { }
};

class Task422 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task422(std::shared_ptr<TATensor<double,6>> I729, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> I730)
   : ta0_(I729), ta1_(Gamma240), ta2_(I730) { }
};

class Task423 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task423(std::shared_ptr<TATensor<double,4>> I730, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I730), ta1_(v2) { }
};

class Task424 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task424(std::shared_ptr<TATensor<double,6>> I729, std::shared_ptr<TATensor<double,6>> Gamma24, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I729), ta1_(Gamma24), ta2_(v2) { }
};

class Task425 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task425(std::shared_ptr<TATensor<double,6>> I729, std::shared_ptr<TATensor<double,6>> Gamma244, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I729), ta1_(Gamma244), ta2_(v2) { }
};

class Task426 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task426(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,6>> I732)
   : ta0_(I72), ta1_(Gamma240), ta2_(I732) { }
};

class Task427 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task427(std::shared_ptr<TATensor<double,6>> I732, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I733)
   : ta0_(I732), ta1_(t2), ta2_(I733) { }
};

class Task428 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task428(std::shared_ptr<TATensor<double,4>> I733, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I733), ta1_(v2) { }
};

class Task429 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task429(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma24, std::shared_ptr<TATensor<double,6>> I738)
   : ta0_(I72), ta1_(Gamma24), ta2_(I738) { }
};

class Task430 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task430(std::shared_ptr<TATensor<double,6>> I738, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I738), ta1_(t2), ta2_(v2) { }
};

class Task431 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task431(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma244, std::shared_ptr<TATensor<double,6>> I744)
   : ta0_(I72), ta1_(Gamma244), ta2_(I744) { }
};

class Task432 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task432(std::shared_ptr<TATensor<double,6>> I744, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I744), ta1_(t2), ta2_(v2) { }
};

class Task433 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task433(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I759)
   : ta0_(I72), ta1_(t2), ta2_(I759) { }
};

class Task434 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task434(std::shared_ptr<TATensor<double,6>> I759, std::shared_ptr<TATensor<double,8>> Gamma250, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I759), ta1_(Gamma250), ta2_(v2) { }
};

class Task435 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task435(std::shared_ptr<TATensor<double,6>> I759, std::shared_ptr<TATensor<double,8>> Gamma251, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I759), ta1_(Gamma251), ta2_(v2) { }
};

class Task436 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task436(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I765)
   : ta0_(I72), ta1_(v2), ta2_(I765) { }
};

class Task437 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task437(std::shared_ptr<TATensor<double,4>> I765, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I765), ta1_(Gamma252), ta2_(t2) { }
};

class Task438 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task438(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I768)
   : ta0_(I72), ta1_(v2), ta2_(I768) { }
};

class Task439 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task439(std::shared_ptr<TATensor<double,4>> I768, std::shared_ptr<TATensor<double,6>> Gamma31, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I768), ta1_(Gamma31), ta2_(t2) { }
};

class Task440 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task440(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I807)
   : ta0_(I72), ta1_(t2), ta2_(I807) { }
};

class Task441 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task441(std::shared_ptr<TATensor<double,4>> I807, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I807), ta1_(Gamma240), ta2_(v2) { }
};

class Task442 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task442(std::shared_ptr<TATensor<double,4>> I807, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I807), ta1_(Gamma252), ta2_(v2) { }
};

class Task443 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task443(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I810)
   : ta0_(I72), ta1_(t2), ta2_(I810) { }
};

class Task444 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task444(std::shared_ptr<TATensor<double,4>> I810, std::shared_ptr<TATensor<double,6>> Gamma230, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I810), ta1_(Gamma230), ta2_(v2) { }
};

class Task445 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task445(std::shared_ptr<TATensor<double,4>> I810, std::shared_ptr<TATensor<double,6>> Gamma31, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I810), ta1_(Gamma31), ta2_(v2) { }
};

class Task446 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task446(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,6>> Gamma276, std::shared_ptr<TATensor<double,6>> I837)
   : ta0_(I72), ta1_(Gamma276), ta2_(I837) { }
};

class Task447 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task447(std::shared_ptr<TATensor<double,6>> I837, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I837), ta1_(t2), ta2_(v2) { }
};

class Task448 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task448(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma568, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I72), ta1_(Gamma568), ta2_(t2) { }
};

class Task449 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task449(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma569, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I72), ta1_(Gamma569), ta2_(t2) { }
};


}
}
}
#endif
#endif

