//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks5.h
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

#ifndef __SRC_SMITH_MRCI_TASKS5_H
#define __SRC_SMITH_MRCI_TASKS5_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task200 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task200(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I37)
   : ta0_(I27), ta1_(Gamma12), ta2_(I37) { }
};

class Task201 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task201(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I37), ta1_(t2), ta2_(h1) { }
};

class Task202 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task202(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I37), ta1_(t2), ta2_(h1) { }
};

class Task203 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task203(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I37), ta1_(t2), ta2_(h1) { }
};

class Task204 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task204(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I37), ta1_(t2), ta2_(h1) { }
};

class Task205 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task205(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I37), ta1_(t2), ta2_(h1) { }
};

class Task206 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task206(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I37), ta1_(t2), ta2_(h1) { }
};

class Task207 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task207(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task208 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task208(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task209 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task209(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task210 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task210(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task211 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task211(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task212 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task212(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task213 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task213(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task214 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task214(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task215 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task215(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task216 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task216(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task217 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task217(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I613)
   : ta0_(I37), ta1_(t2), ta2_(I613) { }
};

class Task218 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task218(std::shared_ptr<TATensor<double,4>> I613, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I613), ta1_(v2) { }
};

class Task219 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task219(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I616)
   : ta0_(I37), ta1_(t2), ta2_(I616) { }
};

class Task220 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task220(std::shared_ptr<TATensor<double,4>> I616, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I616), ta1_(v2) { }
};

class Task221 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task221(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I625)
   : ta0_(I37), ta1_(t2), ta2_(I625) { }
};

class Task222 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task222(std::shared_ptr<TATensor<double,4>> I625, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I625), ta1_(v2) { }
};

class Task223 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task223(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I628)
   : ta0_(I37), ta1_(t2), ta2_(I628) { }
};

class Task224 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task224(std::shared_ptr<TATensor<double,4>> I628, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I628), ta1_(v2) { }
};

class Task225 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task225(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task226 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task226(std::shared_ptr<TATensor<double,4>> I37, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I37), ta1_(t2), ta2_(v2) { }
};

class Task227 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task227(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I55)
   : ta0_(I27), ta1_(h1), ta2_(I55) { }
};

class Task228 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task228(std::shared_ptr<TATensor<double,4>> I55, std::shared_ptr<TATensor<double,4>> Gamma18, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I55), ta1_(Gamma18), ta2_(t2) { }
};

class Task229 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task229(std::shared_ptr<TATensor<double,4>> I55, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I55), ta1_(Gamma10), ta2_(t2) { }
};

class Task230 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task230(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I58)
   : ta0_(I27), ta1_(h1), ta2_(I58) { }
};

class Task231 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task231(std::shared_ptr<TATensor<double,4>> I58, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> I59)
   : ta0_(I58), ta1_(Gamma10), ta2_(I59) { }
};

class Task232 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task232(std::shared_ptr<TATensor<double,4>> I59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I59), ta1_(t2) { }
};

class Task233 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task233(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I67)
   : ta0_(I27), ta1_(t2), ta2_(I67) { }
};

class Task234 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task234(std::shared_ptr<TATensor<double,2>> I67, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I67), ta1_(Gamma12), ta2_(h1) { }
};

class Task235 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task235(std::shared_ptr<TATensor<double,2>> I67, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I67), ta1_(Gamma197), ta2_(v2) { }
};

class Task236 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task236(std::shared_ptr<TATensor<double,2>> I67, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I67), ta1_(Gamma10), ta2_(v2) { }
};

class Task237 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task237(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I70)
   : ta0_(I27), ta1_(t2), ta2_(I70) { }
};

class Task238 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task238(std::shared_ptr<TATensor<double,2>> I70, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I70), ta1_(Gamma12), ta2_(h1) { }
};

class Task239 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task239(std::shared_ptr<TATensor<double,2>> I70, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I70), ta1_(Gamma197), ta2_(v2) { }
};

class Task240 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task240(std::shared_ptr<TATensor<double,2>> I70, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I70), ta1_(Gamma10), ta2_(v2) { }
};

class Task241 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task241(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I387)
   : ta0_(I27), ta1_(t2), ta2_(I387) { }
};

class Task242 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task242(std::shared_ptr<TATensor<double,4>> I387, std::shared_ptr<TATensor<double,6>> Gamma126, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I387), ta1_(Gamma126), ta2_(v2) { }
};

class Task243 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task243(std::shared_ptr<TATensor<double,4>> I387, std::shared_ptr<TATensor<double,6>> Gamma88, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I387), ta1_(Gamma88), ta2_(v2) { }
};

class Task244 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task244(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I393)
   : ta0_(I27), ta1_(v2), ta2_(I393) { }
};

class Task245 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task245(std::shared_ptr<TATensor<double,4>> I393, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I393), ta1_(Gamma0), ta2_(t2) { }
};

class Task246 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task246(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I396)
   : ta0_(I27), ta1_(v2), ta2_(I396) { }
};

class Task247 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task247(std::shared_ptr<TATensor<double,4>> I396, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I396), ta1_(Gamma0), ta2_(t2) { }
};

class Task248 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task248(std::shared_ptr<TATensor<double,4>> I27, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I399)
   : ta0_(I27), ta1_(v2), ta2_(I399) { }
};

class Task249 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task249(std::shared_ptr<TATensor<double,4>> I399, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I399), ta1_(Gamma0), ta2_(t2) { }
};


}
}
}
#endif
#endif

