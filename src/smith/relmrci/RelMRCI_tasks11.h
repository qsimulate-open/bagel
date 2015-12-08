//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks11.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS11_H
#define __SRC_SMITH_RelMRCI_TASKS11_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task500 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task500(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I861)
   : ta0_(I108), ta1_(v2), ta2_(I861) { }
};

class Task501 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task501(std::shared_ptr<TATensor<std::complex<double>,2>> I861, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I861), ta1_(Gamma33), ta2_(t2) { }
};

class Task502 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task502(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I864)
   : ta0_(I108), ta1_(v2), ta2_(I864) { }
};

class Task503 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task503(std::shared_ptr<TATensor<std::complex<double>,2>> I864, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I864), ta1_(Gamma33), ta2_(t2) { }
};

class Task504 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task504(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task505 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task505(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task506 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task506(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task507 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task507(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task508 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task508(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task509 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task509(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task510 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task510(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task511 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task511(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I108), ta1_(t2), ta2_(v2) { }
};

class Task512 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task512(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I935)
   : ta0_(I108), ta1_(t2), ta2_(I935) { }
};

class Task513 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task513(std::shared_ptr<TATensor<std::complex<double>,2>> I935, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I935), ta1_(Gamma24), ta2_(v2) { }
};

class Task514 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task514(std::shared_ptr<TATensor<std::complex<double>,2>> I935, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I935), ta1_(Gamma33), ta2_(v2) { }
};

class Task515 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task515(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I938)
   : ta0_(I108), ta1_(t2), ta2_(I938) { }
};

class Task516 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task516(std::shared_ptr<TATensor<std::complex<double>,2>> I938, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I938), ta1_(Gamma24), ta2_(v2) { }
};

class Task517 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task517(std::shared_ptr<TATensor<std::complex<double>,2>> I938, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I938), ta1_(Gamma33), ta2_(v2) { }
};

class Task518 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task518(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I947)
   : ta0_(I108), ta1_(t2), ta2_(I947) { }
};

class Task519 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task519(std::shared_ptr<TATensor<std::complex<double>,4>> I947, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I947), ta1_(Gamma27), ta2_(v2) { }
};

class Task520 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task520(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I950)
   : ta0_(I108), ta1_(t2), ta2_(I950) { }
};

class Task521 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task521(std::shared_ptr<TATensor<std::complex<double>,4>> I950, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I950), ta1_(Gamma27), ta2_(v2) { }
};

class Task522 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task522(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I953)
   : ta0_(I108), ta1_(t2), ta2_(I953) { }
};

class Task523 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task523(std::shared_ptr<TATensor<std::complex<double>,4>> I953, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I954)
   : ta0_(I953), ta1_(Gamma27), ta2_(I954) { }
};

class Task524 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task524(std::shared_ptr<TATensor<std::complex<double>,4>> I954, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I954), ta1_(v2) { }
};

class Task525 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task525(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I956)
   : ta0_(I108), ta1_(t2), ta2_(I956) { }
};

class Task526 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task526(std::shared_ptr<TATensor<std::complex<double>,4>> I956, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I957)
   : ta0_(I956), ta1_(Gamma27), ta2_(I957) { }
};

class Task527 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task527(std::shared_ptr<TATensor<std::complex<double>,4>> I957, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I957), ta1_(v2) { }
};

class Task528 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task528(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I959)
   : ta0_(I108), ta1_(t2), ta2_(I959) { }
};

class Task529 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task529(std::shared_ptr<TATensor<std::complex<double>,4>> I959, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I960)
   : ta0_(I959), ta1_(Gamma27), ta2_(I960) { }
};

class Task530 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task530(std::shared_ptr<TATensor<std::complex<double>,4>> I960, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I960), ta1_(v2) { }
};

class Task531 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task531(std::shared_ptr<TATensor<std::complex<double>,4>> I108, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I962)
   : ta0_(I108), ta1_(t2), ta2_(I962) { }
};

class Task532 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task532(std::shared_ptr<TATensor<std::complex<double>,4>> I962, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I963)
   : ta0_(I962), ta1_(Gamma27), ta2_(I963) { }
};

class Task533 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task533(std::shared_ptr<TATensor<std::complex<double>,4>> I963, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I963), ta1_(v2) { }
};

class Task534 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task534(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I134)
   : ta0_(proj), ta1_(I134) { }
};

class Task535 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task535(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I135)
   : ta0_(I134), ta1_(h1), ta2_(I135) { }
};

class Task536 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task536(std::shared_ptr<TATensor<std::complex<double>,4>> I135, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I135), ta1_(Gamma24), ta2_(t2) { }
};

class Task537 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task537(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I138)
   : ta0_(I134), ta1_(h1), ta2_(I138) { }
};

class Task538 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task538(std::shared_ptr<TATensor<std::complex<double>,4>> I138, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I138), ta1_(Gamma24), ta2_(t2) { }
};

class Task539 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task539(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,2>> I141)
   : ta0_(I134), ta1_(h1), ta2_(I141) { }
};

class Task540 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task540(std::shared_ptr<TATensor<std::complex<double>,2>> I141, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I141), ta1_(Gamma33), ta2_(t2) { }
};

class Task541 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task541(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,2>> I144)
   : ta0_(I134), ta1_(h1), ta2_(I144) { }
};

class Task542 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task542(std::shared_ptr<TATensor<std::complex<double>,2>> I144, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I144), ta1_(Gamma33), ta2_(t2) { }
};

class Task543 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task543(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I147)
   : ta0_(I134), ta1_(t2), ta2_(I147) { }
};

class Task544 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task544(std::shared_ptr<TATensor<std::complex<double>,2>> I147, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I147), ta1_(Gamma27), ta2_(h1) { }
};

class Task545 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task545(std::shared_ptr<TATensor<std::complex<double>,2>> I147, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I147), ta1_(Gamma33), ta2_(v2) { }
};

class Task546 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task546(std::shared_ptr<TATensor<std::complex<double>,2>> I147, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I147), ta1_(Gamma24), ta2_(v2) { }
};

class Task547 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task547(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> I150)
   : ta0_(I134), ta1_(t2), ta2_(I150) { }
};

class Task548 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task548(std::shared_ptr<TATensor<std::complex<double>,2>> I150, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I150), ta1_(Gamma27), ta2_(h1) { }
};

class Task549 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task549(std::shared_ptr<TATensor<std::complex<double>,2>> I150, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I150), ta1_(Gamma33), ta2_(v2) { }
};


}
}
}
#endif
#endif

