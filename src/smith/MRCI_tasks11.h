//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks11.h
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

#ifndef __SRC_SMITH_MRCI_TASKS11_H
#define __SRC_SMITH_MRCI_TASKS11_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task500 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task500(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I843)
   : ta0_(I108), ta1_(v2), ta2_(I843) { }
};

class Task501 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task501(std::shared_ptr<TATensor<double,6>> I843, std::shared_ptr<TATensor<double,8>> Gamma278, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I843), ta1_(Gamma278), ta2_(t2) { }
};

class Task502 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task502(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I846)
   : ta0_(I108), ta1_(v2), ta2_(I846) { }
};

class Task503 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task503(std::shared_ptr<TATensor<double,6>> I846, std::shared_ptr<TATensor<double,8>> Gamma100, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I846), ta1_(Gamma100), ta2_(t2) { }
};

class Task504 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task504(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I849)
   : ta0_(I108), ta1_(v2), ta2_(I849) { }
};

class Task505 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task505(std::shared_ptr<TATensor<double,4>> I849, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I849), ta1_(Gamma4), ta2_(t2) { }
};

class Task506 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task506(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I852)
   : ta0_(I108), ta1_(v2), ta2_(I852) { }
};

class Task507 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task507(std::shared_ptr<TATensor<double,4>> I852, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I852), ta1_(Gamma4), ta2_(t2) { }
};

class Task508 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task508(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I855)
   : ta0_(I108), ta1_(t2), ta2_(I855) { }
};

class Task509 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task509(std::shared_ptr<TATensor<double,4>> I855, std::shared_ptr<TATensor<double,6>> Gamma221, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I855), ta1_(Gamma221), ta2_(v2) { }
};

class Task510 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task510(std::shared_ptr<TATensor<double,4>> I855, std::shared_ptr<TATensor<double,6>> Gamma104, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I855), ta1_(Gamma104), ta2_(v2) { }
};

class Task511 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task511(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I858)
   : ta0_(I108), ta1_(t2), ta2_(I858) { }
};

class Task512 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task512(std::shared_ptr<TATensor<double,4>> I858, std::shared_ptr<TATensor<double,6>> Gamma221, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I858), ta1_(Gamma221), ta2_(v2) { }
};

class Task513 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task513(std::shared_ptr<TATensor<double,4>> I858, std::shared_ptr<TATensor<double,6>> Gamma104, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I858), ta1_(Gamma104), ta2_(v2) { }
};

class Task514 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task514(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I885)
   : ta0_(I108), ta1_(t2), ta2_(I885) { }
};

class Task515 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task515(std::shared_ptr<TATensor<double,6>> I885, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> I886)
   : ta0_(I885), ta1_(Gamma240), ta2_(I886) { }
};

class Task516 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task516(std::shared_ptr<TATensor<double,4>> I886, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I886), ta1_(v2) { }
};

class Task517 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task517(std::shared_ptr<TATensor<double,6>> I885, std::shared_ptr<TATensor<double,6>> Gamma7, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I885), ta1_(Gamma7), ta2_(v2) { }
};

class Task518 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task518(std::shared_ptr<TATensor<double,6>> I885, std::shared_ptr<TATensor<double,6>> Gamma296, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I885), ta1_(Gamma296), ta2_(v2) { }
};

class Task519 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task519(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,6>> I888)
   : ta0_(I108), ta1_(Gamma240), ta2_(I888) { }
};

class Task520 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task520(std::shared_ptr<TATensor<double,6>> I888, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I889)
   : ta0_(I888), ta1_(t2), ta2_(I889) { }
};

class Task521 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task521(std::shared_ptr<TATensor<double,4>> I889, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I889), ta1_(v2) { }
};

class Task522 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task522(std::shared_ptr<TATensor<double,6>> I888, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I919)
   : ta0_(I888), ta1_(t2), ta2_(I919) { }
};

class Task523 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task523(std::shared_ptr<TATensor<double,4>> I919, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I919), ta1_(v2) { }
};

class Task524 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task524(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,6>> Gamma7, std::shared_ptr<TATensor<double,6>> I894)
   : ta0_(I108), ta1_(Gamma7), ta2_(I894) { }
};

class Task525 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task525(std::shared_ptr<TATensor<double,6>> I894, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I894), ta1_(t2), ta2_(v2) { }
};

class Task526 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task526(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,6>> Gamma296, std::shared_ptr<TATensor<double,6>> I900)
   : ta0_(I108), ta1_(Gamma296), ta2_(I900) { }
};

class Task527 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task527(std::shared_ptr<TATensor<double,6>> I900, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I900), ta1_(t2), ta2_(v2) { }
};

class Task528 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task528(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I915)
   : ta0_(I108), ta1_(t2), ta2_(I915) { }
};

class Task529 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task529(std::shared_ptr<TATensor<double,6>> I915, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> I916)
   : ta0_(I915), ta1_(Gamma240), ta2_(I916) { }
};

class Task530 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task530(std::shared_ptr<TATensor<double,4>> I916, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I916), ta1_(v2) { }
};

class Task531 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task531(std::shared_ptr<TATensor<double,6>> I915, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I915), ta1_(Gamma4), ta2_(v2) { }
};

class Task532 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task532(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,6>> I924)
   : ta0_(I108), ta1_(Gamma4), ta2_(I924) { }
};

class Task533 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task533(std::shared_ptr<TATensor<double,6>> I924, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I924), ta1_(t2), ta2_(v2) { }
};

class Task534 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task534(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I945)
   : ta0_(I108), ta1_(t2), ta2_(I945) { }
};

class Task535 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task535(std::shared_ptr<TATensor<double,6>> I945, std::shared_ptr<TATensor<double,8>> Gamma312, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I945), ta1_(Gamma312), ta2_(v2) { }
};

class Task536 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task536(std::shared_ptr<TATensor<double,6>> I945, std::shared_ptr<TATensor<double,8>> Gamma313, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I945), ta1_(Gamma313), ta2_(v2) { }
};

class Task537 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task537(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I951)
   : ta0_(I108), ta1_(v2), ta2_(I951) { }
};

class Task538 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task538(std::shared_ptr<TATensor<double,4>> I951, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I951), ta1_(Gamma252), ta2_(t2) { }
};

class Task539 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task539(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I954)
   : ta0_(I108), ta1_(v2), ta2_(I954) { }
};

class Task540 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task540(std::shared_ptr<TATensor<double,4>> I954, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I954), ta1_(Gamma252), ta2_(t2) { }
};

class Task541 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task541(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I993)
   : ta0_(I108), ta1_(t2), ta2_(I993) { }
};

class Task542 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task542(std::shared_ptr<TATensor<double,4>> I993, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I993), ta1_(Gamma240), ta2_(v2) { }
};

class Task543 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task543(std::shared_ptr<TATensor<double,4>> I993, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I993), ta1_(Gamma252), ta2_(v2) { }
};

class Task544 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task544(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I996)
   : ta0_(I108), ta1_(t2), ta2_(I996) { }
};

class Task545 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task545(std::shared_ptr<TATensor<double,4>> I996, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I996), ta1_(Gamma240), ta2_(v2) { }
};

class Task546 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task546(std::shared_ptr<TATensor<double,4>> I996, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I996), ta1_(Gamma252), ta2_(v2) { }
};

class Task547 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task547(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,6>> Gamma338, std::shared_ptr<TATensor<double,6>> I1023)
   : ta0_(I108), ta1_(Gamma338), ta2_(I1023) { }
};

class Task548 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task548(std::shared_ptr<TATensor<double,6>> I1023, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1023), ta1_(t2), ta2_(v2) { }
};

class Task549 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task549(std::shared_ptr<TATensor<double,4>> I108, std::shared_ptr<TATensor<double,4>> Gamma569, std::shared_ptr<TATensor<double,4>> I1721)
   : ta0_(I108), ta1_(Gamma569), ta2_(I1721) { }
};


}
}
}
#endif
#endif

