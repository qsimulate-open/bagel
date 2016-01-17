//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks13.h
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

#ifndef __SRC_SMITH_MRCI_TASKS13_H
#define __SRC_SMITH_MRCI_TASKS13_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task600 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task600(std::shared_ptr<TATensor<double,6>> I1107, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1107), ta1_(t2), ta2_(v2) { }
};

class Task601 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task601(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma553, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I144), ta1_(Gamma553), ta2_(t2) { }
};

class Task602 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task602(std::shared_ptr<TATensor<double,4>> I144, std::shared_ptr<TATensor<double,6>> Gamma554, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I144), ta1_(Gamma554), ta2_(t2) { }
};

class Task603 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task603(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I162)
   : ta0_(proj), ta1_(I162) { }
};

class Task604 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task604(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I163)
   : ta0_(I162), ta1_(t2), ta2_(I163) { }
};

class Task605 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task605(std::shared_ptr<TATensor<double,2>> I163, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I163), ta1_(Gamma12), ta2_(h1) { }
};

class Task606 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task606(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I166)
   : ta0_(I162), ta1_(t2), ta2_(I166) { }
};

class Task607 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task607(std::shared_ptr<TATensor<double,2>> I166, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I166), ta1_(Gamma12), ta2_(h1) { }
};

class Task608 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task608(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I169)
   : ta0_(I162), ta1_(h1), ta2_(I169) { }
};

class Task609 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task609(std::shared_ptr<TATensor<double,2>> I169, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I170)
   : ta0_(I169), ta1_(Gamma32), ta2_(I170) { }
};

class Task610 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task610(std::shared_ptr<TATensor<double,4>> I170, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I170), ta1_(t2) { }
};

class Task611 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task611(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,2>> I172)
   : ta0_(I162), ta1_(h1), ta2_(I172) { }
};

class Task612 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task612(std::shared_ptr<TATensor<double,2>> I172, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I173)
   : ta0_(I172), ta1_(Gamma32), ta2_(I173) { }
};

class Task613 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task613(std::shared_ptr<TATensor<double,4>> I173, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I173), ta1_(t2) { }
};

class Task614 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task614(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I181)
   : ta0_(I162), ta1_(t2), ta2_(I181) { }
};

class Task615 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task615(std::shared_ptr<TATensor<double,2>> I181, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I181), ta1_(h1) { }
};

class Task616 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task616(std::shared_ptr<TATensor<double,2>> I181, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1243)
   : ta0_(I181), ta1_(Gamma32), ta2_(I1243) { }
};

class Task617 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task617(std::shared_ptr<TATensor<double,4>> I1243, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1243), ta1_(v2) { }
};

class Task618 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task618(std::shared_ptr<TATensor<double,2>> I181, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I181), ta1_(Gamma12), ta2_(v2) { }
};

class Task619 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task619(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I183)
   : ta0_(I162), ta1_(t2), ta2_(I183) { }
};

class Task620 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task620(std::shared_ptr<TATensor<double,2>> I183, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I183), ta1_(h1) { }
};

class Task621 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task621(std::shared_ptr<TATensor<double,2>> I183, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1246)
   : ta0_(I183), ta1_(Gamma32), ta2_(I1246) { }
};

class Task622 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task622(std::shared_ptr<TATensor<double,4>> I1246, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1246), ta1_(v2) { }
};

class Task623 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task623(std::shared_ptr<TATensor<double,2>> I183, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I183), ta1_(Gamma12), ta2_(v2) { }
};

class Task624 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task624(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I185)
   : ta0_(I162), ta1_(t2), ta2_(I185) { }
};

class Task625 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task625(std::shared_ptr<TATensor<double,2>> I185, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I185), ta1_(h1) { }
};

class Task626 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task626(std::shared_ptr<TATensor<double,2>> I185, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1249)
   : ta0_(I185), ta1_(Gamma32), ta2_(I1249) { }
};

class Task627 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task627(std::shared_ptr<TATensor<double,4>> I1249, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1249), ta1_(v2) { }
};

class Task628 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task628(std::shared_ptr<TATensor<double,2>> I185, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I185), ta1_(Gamma12), ta2_(v2) { }
};

class Task629 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task629(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I187)
   : ta0_(I162), ta1_(t2), ta2_(I187) { }
};

class Task630 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task630(std::shared_ptr<TATensor<double,2>> I187, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I187), ta1_(h1) { }
};

class Task631 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task631(std::shared_ptr<TATensor<double,2>> I187, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1252)
   : ta0_(I187), ta1_(Gamma32), ta2_(I1252) { }
};

class Task632 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task632(std::shared_ptr<TATensor<double,4>> I1252, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1252), ta1_(v2) { }
};

class Task633 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task633(std::shared_ptr<TATensor<double,2>> I187, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I187), ta1_(Gamma12), ta2_(v2) { }
};

class Task634 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task634(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I189)
   : ta0_(I162), ta1_(t2), ta2_(I189) { }
};

class Task635 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task635(std::shared_ptr<TATensor<double,2>> I189, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I189), ta1_(Gamma32), ta2_(h1) { }
};

class Task636 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task636(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I192)
   : ta0_(I162), ta1_(t2), ta2_(I192) { }
};

class Task637 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task637(std::shared_ptr<TATensor<double,2>> I192, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I192), ta1_(Gamma32), ta2_(h1) { }
};

class Task638 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task638(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1116)
   : ta0_(I162), ta1_(v2), ta2_(I1116) { }
};

class Task639 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task639(std::shared_ptr<TATensor<double,2>> I1116, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1116), ta1_(Gamma10), ta2_(t2) { }
};

class Task640 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task640(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1119)
   : ta0_(I162), ta1_(v2), ta2_(I1119) { }
};

class Task641 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task641(std::shared_ptr<TATensor<double,2>> I1119, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1119), ta1_(Gamma10), ta2_(t2) { }
};

class Task642 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task642(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1122)
   : ta0_(I162), ta1_(t2), ta2_(I1122) { }
};

class Task643 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task643(std::shared_ptr<TATensor<double,2>> I1122, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1122), ta1_(Gamma5), ta2_(v2) { }
};

class Task644 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task644(std::shared_ptr<TATensor<double,2>> I1122, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1122), ta1_(Gamma197), ta2_(v2) { }
};

class Task645 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task645(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1125)
   : ta0_(I162), ta1_(t2), ta2_(I1125) { }
};

class Task646 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task646(std::shared_ptr<TATensor<double,2>> I1125, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1125), ta1_(Gamma5), ta2_(v2) { }
};

class Task647 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task647(std::shared_ptr<TATensor<double,2>> I1125, std::shared_ptr<TATensor<double,4>> Gamma197, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1125), ta1_(Gamma197), ta2_(v2) { }
};

class Task648 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task648(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1134)
   : ta0_(I162), ta1_(t2), ta2_(I1134) { }
};

class Task649 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task649(std::shared_ptr<TATensor<double,4>> I1134, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I1135)
   : ta0_(I1134), ta1_(Gamma12), ta2_(I1135) { }
};


}
}
}
#endif
#endif

