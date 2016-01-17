//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks18.h
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

#ifndef __SRC_SMITH_MRCI_TASKS18_H
#define __SRC_SMITH_MRCI_TASKS18_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task850 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task850(std::shared_ptr<TATensor<double,4>> I1511, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1512)
   : ta0_(I1511), ta1_(Gamma51), ta2_(I1512) { }
};

class Task851 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task851(std::shared_ptr<TATensor<double,4>> I1512, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1512), ta1_(v2) { }
};

class Task852 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task852(std::shared_ptr<TATensor<double,4>> I1511, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1511), ta1_(Gamma29), ta2_(v2) { }
};

class Task853 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task853(std::shared_ptr<TATensor<double,4>> I1511, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1511), ta1_(Gamma503), ta2_(v2) { }
};

class Task854 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task854(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1514)
   : ta0_(I194), ta1_(t2), ta2_(I1514) { }
};

class Task855 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task855(std::shared_ptr<TATensor<double,4>> I1514, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1515)
   : ta0_(I1514), ta1_(Gamma51), ta2_(I1515) { }
};

class Task856 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task856(std::shared_ptr<TATensor<double,4>> I1515, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1515), ta1_(v2) { }
};

class Task857 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task857(std::shared_ptr<TATensor<double,4>> I1514, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1514), ta1_(Gamma29), ta2_(v2) { }
};

class Task858 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task858(std::shared_ptr<TATensor<double,4>> I1514, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1514), ta1_(Gamma503), ta2_(v2) { }
};

class Task859 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task859(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1517)
   : ta0_(I194), ta1_(t2), ta2_(I1517) { }
};

class Task860 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task860(std::shared_ptr<TATensor<double,4>> I1517, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1518)
   : ta0_(I1517), ta1_(Gamma51), ta2_(I1518) { }
};

class Task861 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task861(std::shared_ptr<TATensor<double,4>> I1518, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1518), ta1_(v2) { }
};

class Task862 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task862(std::shared_ptr<TATensor<double,4>> I1517, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1517), ta1_(Gamma27), ta2_(v2) { }
};

class Task863 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task863(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1598)
   : ta0_(I194), ta1_(t2), ta2_(I1598) { }
};

class Task864 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task864(std::shared_ptr<TATensor<double,4>> I1598, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1598), ta1_(Gamma50), ta2_(v2) { }
};

class Task865 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task865(std::shared_ptr<TATensor<double,4>> I1598, std::shared_ptr<TATensor<double,6>> Gamma524, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1598), ta1_(Gamma524), ta2_(v2) { }
};

class Task866 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task866(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1604)
   : ta0_(I194), ta1_(v2), ta2_(I1604) { }
};

class Task867 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task867(std::shared_ptr<TATensor<double,4>> I1604, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1604), ta1_(Gamma51), ta2_(t2) { }
};

class Task868 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task868(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1607)
   : ta0_(I194), ta1_(v2), ta2_(I1607) { }
};

class Task869 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task869(std::shared_ptr<TATensor<double,4>> I1607, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1607), ta1_(Gamma51), ta2_(t2) { }
};

class Task870 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task870(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1610)
   : ta0_(I194), ta1_(v2), ta2_(I1610) { }
};

class Task871 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task871(std::shared_ptr<TATensor<double,4>> I1610, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1610), ta1_(Gamma503), ta2_(t2) { }
};

class Task872 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task872(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1613)
   : ta0_(I194), ta1_(v2), ta2_(I1613) { }
};

class Task873 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task873(std::shared_ptr<TATensor<double,4>> I1613, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1613), ta1_(Gamma51), ta2_(t2) { }
};

class Task874 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task874(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> Gamma559, std::shared_ptr<TATensor<double,4>> I1692)
   : ta0_(I194), ta1_(Gamma559), ta2_(I1692) { }
};

class Task875 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task875(std::shared_ptr<TATensor<double,4>> I1692, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1692), ta1_(t2) { }
};

class Task876 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task876(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> Gamma561, std::shared_ptr<TATensor<double,4>> I1696)
   : ta0_(I194), ta1_(Gamma561), ta2_(I1696) { }
};

class Task877 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task877(std::shared_ptr<TATensor<double,4>> I1696, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1696), ta1_(t2) { }
};

class Task878 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task878(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I239)
   : ta0_(proj), ta1_(I239) { }
};

class Task879 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task879(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I240)
   : ta0_(I239), ta1_(h1), ta2_(I240) { }
};

class Task880 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task880(std::shared_ptr<TATensor<double,4>> I240, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I240), ta1_(Gamma50), ta2_(t2) { }
};

class Task881 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task881(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I243)
   : ta0_(I239), ta1_(Gamma51), ta2_(I243) { }
};

class Task882 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task882(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I243), ta1_(t2), ta2_(h1) { }
};

class Task883 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task883(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I243), ta1_(t2), ta2_(h1) { }
};

class Task884 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task884(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task885 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task885(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task886 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task886(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task887 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task887(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1616)
   : ta0_(I239), ta1_(t2), ta2_(I1616) { }
};

class Task888 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task888(std::shared_ptr<TATensor<double,6>> I1616, std::shared_ptr<TATensor<double,6>> Gamma529, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1616), ta1_(Gamma529), ta2_(v2) { }
};

class Task889 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task889(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1619)
   : ta0_(I239), ta1_(t2), ta2_(I1619) { }
};

class Task890 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task890(std::shared_ptr<TATensor<double,6>> I1619, std::shared_ptr<TATensor<double,6>> Gamma530, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1619), ta1_(Gamma530), ta2_(v2) { }
};

class Task891 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task891(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1622)
   : ta0_(I239), ta1_(t2), ta2_(I1622) { }
};

class Task892 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task892(std::shared_ptr<TATensor<double,6>> I1622, std::shared_ptr<TATensor<double,8>> Gamma531, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1622), ta1_(Gamma531), ta2_(v2) { }
};

class Task893 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task893(std::shared_ptr<TATensor<double,6>> I1622, std::shared_ptr<TATensor<double,8>> Gamma349, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1622), ta1_(Gamma349), ta2_(v2) { }
};

class Task894 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task894(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1628)
   : ta0_(I239), ta1_(v2), ta2_(I1628) { }
};

class Task895 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task895(std::shared_ptr<TATensor<double,4>> I1628, std::shared_ptr<TATensor<double,6>> Gamma471, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1628), ta1_(Gamma471), ta2_(t2) { }
};

class Task896 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task896(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1634)
   : ta0_(I239), ta1_(t2), ta2_(I1634) { }
};

class Task897 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task897(std::shared_ptr<TATensor<double,4>> I1634, std::shared_ptr<TATensor<double,6>> Gamma524, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1634), ta1_(Gamma524), ta2_(v2) { }
};

class Task898 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task898(std::shared_ptr<TATensor<double,4>> I1634, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1634), ta1_(Gamma50), ta2_(v2) { }
};

class Task899 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task899(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> I1640)
   : ta0_(I239), ta1_(Gamma503), ta2_(I1640) { }
};


}
}
}
#endif
#endif

