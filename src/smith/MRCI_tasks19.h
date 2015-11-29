//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks19.h
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

#ifndef __SRC_SMITH_MRCI_TASKS19_H
#define __SRC_SMITH_MRCI_TASKS19_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task900 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task900(std::shared_ptr<TATensor<double,4>> I1644, std::shared_ptr<TATensor<double,6>> Gamma526, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1644), ta1_(Gamma526), ta2_(v2) { }
};

class Task901 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task901(std::shared_ptr<TATensor<double,4>> I1644, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1644), ta1_(Gamma50), ta2_(v2) { }
};

class Task902 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task902(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> I1650)
   : ta0_(I239), ta1_(Gamma503), ta2_(I1650) { }
};

class Task903 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task903(std::shared_ptr<TATensor<double,4>> I1650, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1650), ta1_(t2), ta2_(v2) { }
};

class Task904 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task904(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,6>> Gamma526, std::shared_ptr<TATensor<double,6>> I1662)
   : ta0_(I239), ta1_(Gamma526), ta2_(I1662) { }
};

class Task905 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task905(std::shared_ptr<TATensor<double,6>> I1662, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1663)
   : ta0_(I1662), ta1_(t2), ta2_(I1663) { }
};

class Task906 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task906(std::shared_ptr<TATensor<double,4>> I1663, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1663), ta1_(v2) { }
};

class Task907 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task907(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,6>> I1665)
   : ta0_(I239), ta1_(Gamma50), ta2_(I1665) { }
};

class Task908 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task908(std::shared_ptr<TATensor<double,6>> I1665, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1665), ta1_(t2), ta2_(v2) { }
};

class Task909 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task909(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,6>> Gamma545, std::shared_ptr<TATensor<double,6>> I1668)
   : ta0_(I239), ta1_(Gamma545), ta2_(I1668) { }
};

class Task910 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task910(std::shared_ptr<TATensor<double,6>> I1668, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1668), ta1_(t2), ta2_(v2) { }
};

class Task911 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task911(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I260)
   : ta0_(proj), ta1_(I260) { }
};

class Task912 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task912(std::shared_ptr<TATensor<double,4>> I260, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> I261)
   : ta0_(I260), ta1_(Gamma2), ta2_(I261) { }
};

class Task913 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task913(std::shared_ptr<TATensor<double,4>> I261, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I261), ta1_(t2), ta2_(v2) { }
};

class Task914 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task914(std::shared_ptr<TATensor<double,4>> I261, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I261), ta1_(t2), ta2_(v2) { }
};

class Task915 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task915(std::shared_ptr<TATensor<double,4>> I260, std::shared_ptr<TATensor<double,4>> Gamma548, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I260), ta1_(Gamma548), ta2_(t2) { }
};

class Task916 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task916(std::shared_ptr<TATensor<double,4>> I260, std::shared_ptr<TATensor<double,4>> Gamma549, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I260), ta1_(Gamma549), ta2_(t2) { }
};

class Task917 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task917(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1112)
   : ta0_(proj), ta1_(I1112) { }
};

class Task918 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task918(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1113)
   : ta0_(I1112), ta1_(v2), ta2_(I1113) { }
};

class Task919 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task919(std::shared_ptr<TATensor<double,4>> I1113, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1113), ta1_(Gamma2), ta2_(t2) { }
};

class Task920 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task920(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1112), ta1_(t2), ta2_(v2) { }
};

class Task921 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task921(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1112), ta1_(t2), ta2_(v2) { }
};

class Task922 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task922(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1112), ta1_(t2), ta2_(v2) { }
};

class Task923 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task923(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1112), ta1_(t2), ta2_(v2) { }
};

class Task924 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task924(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1356)
   : ta0_(I1112), ta1_(t2), ta2_(I1356) { }
};

class Task925 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task925(std::shared_ptr<TATensor<double,4>> I1356, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1356), ta1_(Gamma503), ta2_(v2) { }
};

class Task926 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,0>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task926(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,0>> Gamma558, std::shared_ptr<TATensor<double,4>> I1697)
   : ta0_(I1112), ta1_(Gamma558), ta2_(I1697) { }
};

class Task927 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task927(std::shared_ptr<TATensor<double,4>> I1697, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1697), ta1_(t2) { }
};

class Task928 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,0>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task928(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,0>> Gamma560, std::shared_ptr<TATensor<double,4>> I1701)
   : ta0_(I1112), ta1_(Gamma560), ta2_(I1701) { }
};

class Task929 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task929(std::shared_ptr<TATensor<double,4>> I1701, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1701), ta1_(t2) { }
};

class Task930 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task930(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1640)
   : ta0_(proj), ta1_(I1640) { }
};

class Task931 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task931(std::shared_ptr<TATensor<double,4>> I1640, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1641)
   : ta0_(I1640), ta1_(t2), ta2_(I1641) { }
};

class Task932 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task932(std::shared_ptr<TATensor<double,4>> I1641, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1641), ta1_(Gamma503), ta2_(v2) { }
};

class Task933 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task933(std::shared_ptr<TATensor<double,4>> I1640, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> I1674)
   : ta0_(I1640), ta1_(Gamma503), ta2_(I1674) { }
};

class Task934 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task934(std::shared_ptr<TATensor<double,4>> I1674, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1674), ta1_(t2), ta2_(v2) { }
};

class Task935 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task935(std::shared_ptr<TATensor<double,4>> I1640, std::shared_ptr<TATensor<double,4>> Gamma566, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1640), ta1_(Gamma566), ta2_(t2) { }
};

class Task936 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task936(std::shared_ptr<TATensor<double,4>> I1640, std::shared_ptr<TATensor<double,4>> Gamma567, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1640), ta1_(Gamma567), ta2_(t2) { }
};

class Task937 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task937(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task938 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task938(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1732)
   : ta0_(proj), ta1_(I1732) { }
};

class Task939 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task939(std::shared_ptr<TATensor<double,4>> I1732, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1732), ta1_(Gamma5), ta2_(h1) { }
};

class Task940 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task940(std::shared_ptr<TATensor<double,4>> I1732, std::shared_ptr<TATensor<double,6>> Gamma104, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1732), ta1_(Gamma104), ta2_(v2) { }
};

class Task941 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task941(std::shared_ptr<TATensor<double,4>> I1732, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1732), ta1_(Gamma4), ta2_(v2) { }
};

class Task942 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task942(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1734)
   : ta0_(proj), ta1_(I1734) { }
};

class Task943 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task943(std::shared_ptr<TATensor<double,4>> I1734, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1734), ta1_(Gamma32), ta2_(h1) { }
};

class Task944 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task944(std::shared_ptr<TATensor<double,4>> I1734, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1751)
   : ta0_(I1734), ta1_(Gamma29), ta2_(I1751) { }
};

class Task945 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task945(std::shared_ptr<TATensor<double,4>> I1751, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1751), ta1_(v2) { }
};

class Task946 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task946(std::shared_ptr<TATensor<double,4>> I1734, std::shared_ptr<TATensor<double,4>> Gamma25, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1734), ta1_(Gamma25), ta2_(v2) { }
};

class Task947 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task947(std::shared_ptr<TATensor<double,4>> I1734, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1734), ta1_(Gamma27), ta2_(v2) { }
};

class Task948 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task948(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1736)
   : ta0_(proj), ta1_(I1736) { }
};

class Task949 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task949(std::shared_ptr<TATensor<double,4>> I1736, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1736), ta1_(Gamma32), ta2_(h1) { }
};


}
}
}
#endif
#endif

