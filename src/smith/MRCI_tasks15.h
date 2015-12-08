//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks15.h
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

#ifndef __SRC_SMITH_MRCI_TASKS15_H
#define __SRC_SMITH_MRCI_TASKS15_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task700 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task700(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task701 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task701(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task702 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task702(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task703 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task703(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task704 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task704(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task705 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task705(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1310)
   : ta0_(I162), ta1_(t2), ta2_(I1310) { }
};

class Task706 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task706(std::shared_ptr<TATensor<double,2>> I1310, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1310), ta1_(Gamma29), ta2_(v2) { }
};

class Task707 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task707(std::shared_ptr<TATensor<double,2>> I1310, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1310), ta1_(Gamma51), ta2_(v2) { }
};

class Task708 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task708(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1313)
   : ta0_(I162), ta1_(t2), ta2_(I1313) { }
};

class Task709 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task709(std::shared_ptr<TATensor<double,2>> I1313, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1313), ta1_(Gamma29), ta2_(v2) { }
};

class Task710 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task710(std::shared_ptr<TATensor<double,2>> I1313, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1313), ta1_(Gamma51), ta2_(v2) { }
};

class Task711 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task711(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1322)
   : ta0_(I162), ta1_(t2), ta2_(I1322) { }
};

class Task712 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task712(std::shared_ptr<TATensor<double,4>> I1322, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1322), ta1_(Gamma32), ta2_(v2) { }
};

class Task713 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task713(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1325)
   : ta0_(I162), ta1_(t2), ta2_(I1325) { }
};

class Task714 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task714(std::shared_ptr<TATensor<double,4>> I1325, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1325), ta1_(Gamma32), ta2_(v2) { }
};

class Task715 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task715(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1328)
   : ta0_(I162), ta1_(t2), ta2_(I1328) { }
};

class Task716 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task716(std::shared_ptr<TATensor<double,4>> I1328, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1329)
   : ta0_(I1328), ta1_(Gamma32), ta2_(I1329) { }
};

class Task717 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task717(std::shared_ptr<TATensor<double,4>> I1329, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1329), ta1_(v2) { }
};

class Task718 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task718(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1331)
   : ta0_(I162), ta1_(t2), ta2_(I1331) { }
};

class Task719 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task719(std::shared_ptr<TATensor<double,4>> I1331, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1332)
   : ta0_(I1331), ta1_(Gamma32), ta2_(I1332) { }
};

class Task720 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task720(std::shared_ptr<TATensor<double,4>> I1332, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1332), ta1_(v2) { }
};

class Task721 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task721(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1334)
   : ta0_(I162), ta1_(t2), ta2_(I1334) { }
};

class Task722 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task722(std::shared_ptr<TATensor<double,4>> I1334, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1335)
   : ta0_(I1334), ta1_(Gamma32), ta2_(I1335) { }
};

class Task723 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task723(std::shared_ptr<TATensor<double,4>> I1335, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1335), ta1_(v2) { }
};

class Task724 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task724(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1337)
   : ta0_(I162), ta1_(t2), ta2_(I1337) { }
};

class Task725 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task725(std::shared_ptr<TATensor<double,4>> I1337, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1338)
   : ta0_(I1337), ta1_(Gamma32), ta2_(I1338) { }
};

class Task726 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task726(std::shared_ptr<TATensor<double,4>> I1338, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1338), ta1_(v2) { }
};

class Task727 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task727(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I194)
   : ta0_(proj), ta1_(I194) { }
};

class Task728 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task728(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I195)
   : ta0_(I194), ta1_(h1), ta2_(I195) { }
};

class Task729 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task729(std::shared_ptr<TATensor<double,4>> I195, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I196)
   : ta0_(I195), ta1_(Gamma29), ta2_(I196) { }
};

class Task730 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task730(std::shared_ptr<TATensor<double,4>> I196, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I196), ta1_(t2) { }
};

class Task731 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task731(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I198)
   : ta0_(I194), ta1_(h1), ta2_(I198) { }
};

class Task732 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task732(std::shared_ptr<TATensor<double,4>> I198, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I198), ta1_(Gamma27), ta2_(t2) { }
};

class Task733 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task733(std::shared_ptr<TATensor<double,4>> I198, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I198), ta1_(Gamma29), ta2_(t2) { }
};

class Task734 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task734(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I207)
   : ta0_(I194), ta1_(h1), ta2_(I207) { }
};

class Task735 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task735(std::shared_ptr<TATensor<double,2>> I207, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I207), ta1_(Gamma51), ta2_(t2) { }
};

class Task736 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task736(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I210)
   : ta0_(I194), ta1_(h1), ta2_(I210) { }
};

class Task737 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task737(std::shared_ptr<TATensor<double,2>> I210, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I210), ta1_(Gamma51), ta2_(t2) { }
};

class Task738 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task738(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I213)
   : ta0_(I194), ta1_(t2), ta2_(I213) { }
};

class Task739 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task739(std::shared_ptr<TATensor<double,2>> I213, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I213), ta1_(Gamma32), ta2_(h1) { }
};

class Task740 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task740(std::shared_ptr<TATensor<double,2>> I213, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I213), ta1_(Gamma51), ta2_(v2) { }
};

class Task741 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task741(std::shared_ptr<TATensor<double,2>> I213, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I213), ta1_(Gamma29), ta2_(v2) { }
};

class Task742 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task742(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I216)
   : ta0_(I194), ta1_(t2), ta2_(I216) { }
};

class Task743 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task743(std::shared_ptr<TATensor<double,2>> I216, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I216), ta1_(Gamma32), ta2_(h1) { }
};

class Task744 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task744(std::shared_ptr<TATensor<double,2>> I216, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I216), ta1_(Gamma51), ta2_(v2) { }
};

class Task745 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task745(std::shared_ptr<TATensor<double,2>> I216, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I216), ta1_(Gamma29), ta2_(v2) { }
};

class Task746 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task746(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I219)
   : ta0_(I194), ta1_(Gamma32), ta2_(I219) { }
};

class Task747 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task747(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I219), ta1_(t2), ta2_(h1) { }
};

class Task748 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task748(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I219), ta1_(t2), ta2_(h1) { }
};

class Task749 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task749(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I219), ta1_(t2), ta2_(h1) { }
};


}
}
}
#endif
#endif

