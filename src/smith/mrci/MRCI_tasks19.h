//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks19.h
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
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task900(std::shared_ptr<TATensor<double,4>> I1640, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1640), ta1_(t2), ta2_(v2) { }
};

class Task901 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task901(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,6>> Gamma524, std::shared_ptr<TATensor<double,6>> I1652)
   : ta0_(I239), ta1_(Gamma524), ta2_(I1652) { }
};

class Task902 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task902(std::shared_ptr<TATensor<double,6>> I1652, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1653)
   : ta0_(I1652), ta1_(t2), ta2_(I1653) { }
};

class Task903 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task903(std::shared_ptr<TATensor<double,4>> I1653, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1653), ta1_(v2) { }
};

class Task904 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task904(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,6>> I1655)
   : ta0_(I239), ta1_(Gamma50), ta2_(I1655) { }
};

class Task905 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task905(std::shared_ptr<TATensor<double,6>> I1655, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1655), ta1_(t2), ta2_(v2) { }
};

class Task906 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task906(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,6>> Gamma543, std::shared_ptr<TATensor<double,6>> I1658)
   : ta0_(I239), ta1_(Gamma543), ta2_(I1658) { }
};

class Task907 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task907(std::shared_ptr<TATensor<double,6>> I1658, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1658), ta1_(t2), ta2_(v2) { }
};

class Task908 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task908(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I260)
   : ta0_(proj), ta1_(I260) { }
};

class Task909 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task909(std::shared_ptr<TATensor<double,4>> I260, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> I261)
   : ta0_(I260), ta1_(Gamma2), ta2_(I261) { }
};

class Task910 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task910(std::shared_ptr<TATensor<double,4>> I261, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I261), ta1_(t2), ta2_(v2) { }
};

class Task911 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task911(std::shared_ptr<TATensor<double,4>> I261, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I261), ta1_(t2), ta2_(v2) { }
};

class Task912 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task912(std::shared_ptr<TATensor<double,4>> I260, std::shared_ptr<TATensor<double,4>> Gamma545, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I260), ta1_(Gamma545), ta2_(t2) { }
};

class Task913 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task913(std::shared_ptr<TATensor<double,4>> I260, std::shared_ptr<TATensor<double,4>> Gamma546, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I260), ta1_(Gamma546), ta2_(t2) { }
};

class Task914 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task914(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1112)
   : ta0_(proj), ta1_(I1112) { }
};

class Task915 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task915(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1113)
   : ta0_(I1112), ta1_(v2), ta2_(I1113) { }
};

class Task916 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task916(std::shared_ptr<TATensor<double,4>> I1113, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1113), ta1_(Gamma2), ta2_(t2) { }
};

class Task917 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task917(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1112), ta1_(t2), ta2_(v2) { }
};

class Task918 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task918(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1112), ta1_(t2), ta2_(v2) { }
};

class Task919 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task919(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1352)
   : ta0_(I1112), ta1_(t2), ta2_(I1352) { }
};

class Task920 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task920(std::shared_ptr<TATensor<double,4>> I1352, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1352), ta1_(Gamma503), ta2_(v2) { }
};

class Task921 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,0>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task921(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,0>> Gamma555, std::shared_ptr<TATensor<double,4>> I1684)
   : ta0_(I1112), ta1_(Gamma555), ta2_(I1684) { }
};

class Task922 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task922(std::shared_ptr<TATensor<double,4>> I1684, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1684), ta1_(t2) { }
};

class Task923 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,0>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task923(std::shared_ptr<TATensor<double,4>> I1112, std::shared_ptr<TATensor<double,0>> Gamma557, std::shared_ptr<TATensor<double,4>> I1688)
   : ta0_(I1112), ta1_(Gamma557), ta2_(I1688) { }
};

class Task924 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task924(std::shared_ptr<TATensor<double,4>> I1688, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1688), ta1_(t2) { }
};

class Task925 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task925(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1630)
   : ta0_(proj), ta1_(I1630) { }
};

class Task926 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task926(std::shared_ptr<TATensor<double,4>> I1630, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1631)
   : ta0_(I1630), ta1_(t2), ta2_(I1631) { }
};

class Task927 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task927(std::shared_ptr<TATensor<double,4>> I1631, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1631), ta1_(Gamma503), ta2_(v2) { }
};

class Task928 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task928(std::shared_ptr<TATensor<double,4>> I1630, std::shared_ptr<TATensor<double,4>> Gamma563, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1630), ta1_(Gamma563), ta2_(t2) { }
};

class Task929 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task929(std::shared_ptr<TATensor<double,4>> I1630, std::shared_ptr<TATensor<double,4>> Gamma564, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1630), ta1_(Gamma564), ta2_(t2) { }
};

class Task930 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task930(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task931 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task931(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1719)
   : ta0_(proj), ta1_(I1719) { }
};

class Task932 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task932(std::shared_ptr<TATensor<double,4>> I1719, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1719), ta1_(Gamma5), ta2_(h1) { }
};

class Task933 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task933(std::shared_ptr<TATensor<double,4>> I1719, std::shared_ptr<TATensor<double,6>> Gamma104, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1719), ta1_(Gamma104), ta2_(v2) { }
};

class Task934 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task934(std::shared_ptr<TATensor<double,4>> I1719, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1719), ta1_(Gamma4), ta2_(v2) { }
};

class Task935 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task935(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1721)
   : ta0_(proj), ta1_(I1721) { }
};

class Task936 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task936(std::shared_ptr<TATensor<double,4>> I1721, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1721), ta1_(Gamma32), ta2_(h1) { }
};

class Task937 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task937(std::shared_ptr<TATensor<double,4>> I1721, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1738)
   : ta0_(I1721), ta1_(Gamma29), ta2_(I1738) { }
};

class Task938 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task938(std::shared_ptr<TATensor<double,4>> I1738, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1738), ta1_(v2) { }
};

class Task939 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task939(std::shared_ptr<TATensor<double,4>> I1721, std::shared_ptr<TATensor<double,4>> Gamma25, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1721), ta1_(Gamma25), ta2_(v2) { }
};

class Task940 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task940(std::shared_ptr<TATensor<double,4>> I1721, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1721), ta1_(Gamma27), ta2_(v2) { }
};

class Task941 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task941(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1723)
   : ta0_(proj), ta1_(I1723) { }
};

class Task942 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task942(std::shared_ptr<TATensor<double,4>> I1723, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1723), ta1_(Gamma32), ta2_(h1) { }
};

class Task943 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task943(std::shared_ptr<TATensor<double,4>> I1723, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1746)
   : ta0_(I1723), ta1_(Gamma29), ta2_(I1746) { }
};

class Task944 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task944(std::shared_ptr<TATensor<double,4>> I1746, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1746), ta1_(v2) { }
};

class Task945 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task945(std::shared_ptr<TATensor<double,4>> I1723, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1723), ta1_(Gamma5), ta2_(v2) { }
};

class Task946 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task946(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1725)
   : ta0_(proj), ta1_(I1725) { }
};

class Task947 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task947(std::shared_ptr<TATensor<double,4>> I1725, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1725), ta1_(Gamma51), ta2_(h1) { }
};

class Task948 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task948(std::shared_ptr<TATensor<double,4>> I1725, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1725), ta1_(Gamma50), ta2_(v2) { }
};

class Task949 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task949(std::shared_ptr<TATensor<double,4>> I1725, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1725), ta1_(Gamma49), ta2_(v2) { }
};


}
}
}
#endif
#endif

