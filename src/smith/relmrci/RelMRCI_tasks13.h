//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks13.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS13_H
#define __SRC_SMITH_RelMRCI_TASKS13_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task600 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task600(std::shared_ptr<TATensor<std::complex<double>,4>> I1016, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1016), ta1_(Gamma24), ta2_(t2) { }
};

class Task601 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task601(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1019)
   : ta0_(I134), ta1_(v2), ta2_(I1019) { }
};

class Task602 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task602(std::shared_ptr<TATensor<std::complex<double>,4>> I1019, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1019), ta1_(Gamma24), ta2_(t2) { }
};

class Task603 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task603(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1022)
   : ta0_(I134), ta1_(v2), ta2_(I1022) { }
};

class Task604 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task604(std::shared_ptr<TATensor<std::complex<double>,4>> I1022, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1022), ta1_(Gamma24), ta2_(t2) { }
};

class Task605 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task605(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1025)
   : ta0_(I134), ta1_(v2), ta2_(I1025) { }
};

class Task606 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task606(std::shared_ptr<TATensor<std::complex<double>,4>> I1025, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1025), ta1_(Gamma32), ta2_(t2) { }
};

class Task607 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task607(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1028)
   : ta0_(I134), ta1_(v2), ta2_(I1028) { }
};

class Task608 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task608(std::shared_ptr<TATensor<std::complex<double>,4>> I1028, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1028), ta1_(Gamma32), ta2_(t2) { }
};

class Task609 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task609(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1031)
   : ta0_(I134), ta1_(v2), ta2_(I1031) { }
};

class Task610 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task610(std::shared_ptr<TATensor<std::complex<double>,4>> I1031, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma26, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1031), ta1_(Gamma26), ta2_(t2) { }
};

class Task611 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task611(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1034)
   : ta0_(I134), ta1_(v2), ta2_(I1034) { }
};

class Task612 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task612(std::shared_ptr<TATensor<std::complex<double>,4>> I1034, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma335, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1034), ta1_(Gamma335), ta2_(t2) { }
};

class Task613 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task613(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1037)
   : ta0_(I134), ta1_(v2), ta2_(I1037) { }
};

class Task614 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task614(std::shared_ptr<TATensor<std::complex<double>,4>> I1037, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma336, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1037), ta1_(Gamma336), ta2_(t2) { }
};

class Task615 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task615(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1040)
   : ta0_(I134), ta1_(v2), ta2_(I1040) { }
};

class Task616 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task616(std::shared_ptr<TATensor<std::complex<double>,4>> I1040, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1040), ta1_(Gamma32), ta2_(t2) { }
};

class Task617 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task617(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1043)
   : ta0_(I134), ta1_(v2), ta2_(I1043) { }
};

class Task618 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task618(std::shared_ptr<TATensor<std::complex<double>,4>> I1043, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1043), ta1_(Gamma32), ta2_(t2) { }
};

class Task619 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task619(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1046)
   : ta0_(I134), ta1_(v2), ta2_(I1046) { }
};

class Task620 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task620(std::shared_ptr<TATensor<std::complex<double>,4>> I1046, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1046), ta1_(Gamma32), ta2_(t2) { }
};

class Task621 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task621(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I1049)
   : ta0_(I134), ta1_(v2), ta2_(I1049) { }
};

class Task622 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task622(std::shared_ptr<TATensor<std::complex<double>,2>> I1049, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1049), ta1_(Gamma33), ta2_(t2) { }
};

class Task623 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task623(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,2>> I1052)
   : ta0_(I134), ta1_(v2), ta2_(I1052) { }
};

class Task624 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task624(std::shared_ptr<TATensor<std::complex<double>,2>> I1052, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1052), ta1_(Gamma33), ta2_(t2) { }
};

class Task625 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task625(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1067)
   : ta0_(I134), ta1_(t2), ta2_(I1067) { }
};

class Task626 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task626(std::shared_ptr<TATensor<std::complex<double>,4>> I1067, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1067), ta1_(Gamma27), ta2_(v2) { }
};

class Task627 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task627(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1070)
   : ta0_(I134), ta1_(t2), ta2_(I1070) { }
};

class Task628 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task628(std::shared_ptr<TATensor<std::complex<double>,4>> I1070, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1070), ta1_(Gamma27), ta2_(v2) { }
};

class Task629 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task629(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1097)
   : ta0_(I134), ta1_(t2), ta2_(I1097) { }
};

class Task630 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task630(std::shared_ptr<TATensor<std::complex<double>,4>> I1097, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1098)
   : ta0_(I1097), ta1_(Gamma33), ta2_(I1098) { }
};

class Task631 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task631(std::shared_ptr<TATensor<std::complex<double>,4>> I1098, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1098), ta1_(v2) { }
};

class Task632 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task632(std::shared_ptr<TATensor<std::complex<double>,4>> I1097, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1097), ta1_(Gamma24), ta2_(v2) { }
};

class Task633 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task633(std::shared_ptr<TATensor<std::complex<double>,4>> I1097, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1097), ta1_(Gamma368), ta2_(v2) { }
};

class Task634 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task634(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1100)
   : ta0_(I134), ta1_(t2), ta2_(I1100) { }
};

class Task635 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task635(std::shared_ptr<TATensor<std::complex<double>,4>> I1100, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1101)
   : ta0_(I1100), ta1_(Gamma33), ta2_(I1101) { }
};

class Task636 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task636(std::shared_ptr<TATensor<std::complex<double>,4>> I1101, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1101), ta1_(v2) { }
};

class Task637 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task637(std::shared_ptr<TATensor<std::complex<double>,4>> I1100, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma363, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1100), ta1_(Gamma363), ta2_(v2) { }
};

class Task638 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task638(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1103)
   : ta0_(I134), ta1_(t2), ta2_(I1103) { }
};

class Task639 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task639(std::shared_ptr<TATensor<std::complex<double>,4>> I1103, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1104)
   : ta0_(I1103), ta1_(Gamma33), ta2_(I1104) { }
};

class Task640 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task640(std::shared_ptr<TATensor<std::complex<double>,4>> I1104, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1104), ta1_(v2) { }
};

class Task641 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task641(std::shared_ptr<TATensor<std::complex<double>,4>> I1103, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1103), ta1_(Gamma24), ta2_(v2) { }
};

class Task642 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task642(std::shared_ptr<TATensor<std::complex<double>,4>> I1103, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1103), ta1_(Gamma368), ta2_(v2) { }
};

class Task643 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task643(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1106)
   : ta0_(I134), ta1_(t2), ta2_(I1106) { }
};

class Task644 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task644(std::shared_ptr<TATensor<std::complex<double>,4>> I1106, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1107)
   : ta0_(I1106), ta1_(Gamma33), ta2_(I1107) { }
};

class Task645 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task645(std::shared_ptr<TATensor<std::complex<double>,4>> I1107, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1107), ta1_(v2) { }
};

class Task646 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task646(std::shared_ptr<TATensor<std::complex<double>,4>> I1106, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1106), ta1_(Gamma24), ta2_(v2) { }
};

class Task647 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task647(std::shared_ptr<TATensor<std::complex<double>,4>> I1106, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1106), ta1_(Gamma368), ta2_(v2) { }
};

class Task648 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task648(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1109)
   : ta0_(I134), ta1_(t2), ta2_(I1109) { }
};

class Task649 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task649(std::shared_ptr<TATensor<std::complex<double>,4>> I1109, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1110)
   : ta0_(I1109), ta1_(Gamma33), ta2_(I1110) { }
};


}
}
}
#endif
#endif

