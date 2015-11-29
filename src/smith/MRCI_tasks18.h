//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks18.h
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
    Task850(std::shared_ptr<TATensor<double,4>> I1512, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1512), ta1_(Gamma29), ta2_(v2) { }
};

class Task851 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task851(std::shared_ptr<TATensor<double,4>> I1512, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1512), ta1_(Gamma503), ta2_(v2) { }
};

class Task852 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task852(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1515)
   : ta0_(I194), ta1_(t2), ta2_(I1515) { }
};

class Task853 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task853(std::shared_ptr<TATensor<double,4>> I1515, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1516)
   : ta0_(I1515), ta1_(Gamma51), ta2_(I1516) { }
};

class Task854 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task854(std::shared_ptr<TATensor<double,4>> I1516, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1516), ta1_(v2) { }
};

class Task855 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task855(std::shared_ptr<TATensor<double,4>> I1515, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1515), ta1_(Gamma29), ta2_(v2) { }
};

class Task856 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task856(std::shared_ptr<TATensor<double,4>> I1515, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1515), ta1_(Gamma503), ta2_(v2) { }
};

class Task857 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task857(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1518)
   : ta0_(I194), ta1_(t2), ta2_(I1518) { }
};

class Task858 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task858(std::shared_ptr<TATensor<double,4>> I1518, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1519)
   : ta0_(I1518), ta1_(Gamma51), ta2_(I1519) { }
};

class Task859 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task859(std::shared_ptr<TATensor<double,4>> I1519, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1519), ta1_(v2) { }
};

class Task860 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task860(std::shared_ptr<TATensor<double,4>> I1518, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1518), ta1_(Gamma29), ta2_(v2) { }
};

class Task861 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task861(std::shared_ptr<TATensor<double,4>> I1518, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1518), ta1_(Gamma503), ta2_(v2) { }
};

class Task862 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task862(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1521)
   : ta0_(I194), ta1_(t2), ta2_(I1521) { }
};

class Task863 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task863(std::shared_ptr<TATensor<double,4>> I1521, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1522)
   : ta0_(I1521), ta1_(Gamma51), ta2_(I1522) { }
};

class Task864 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task864(std::shared_ptr<TATensor<double,4>> I1522, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1522), ta1_(v2) { }
};

class Task865 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task865(std::shared_ptr<TATensor<double,4>> I1521, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1521), ta1_(Gamma27), ta2_(v2) { }
};

class Task866 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task866(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1608)
   : ta0_(I194), ta1_(t2), ta2_(I1608) { }
};

class Task867 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task867(std::shared_ptr<TATensor<double,4>> I1608, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1608), ta1_(Gamma50), ta2_(v2) { }
};

class Task868 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task868(std::shared_ptr<TATensor<double,4>> I1608, std::shared_ptr<TATensor<double,6>> Gamma526, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1608), ta1_(Gamma526), ta2_(v2) { }
};

class Task869 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task869(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1614)
   : ta0_(I194), ta1_(v2), ta2_(I1614) { }
};

class Task870 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task870(std::shared_ptr<TATensor<double,4>> I1614, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1614), ta1_(Gamma51), ta2_(t2) { }
};

class Task871 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task871(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1617)
   : ta0_(I194), ta1_(v2), ta2_(I1617) { }
};

class Task872 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task872(std::shared_ptr<TATensor<double,4>> I1617, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1617), ta1_(Gamma51), ta2_(t2) { }
};

class Task873 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task873(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1620)
   : ta0_(I194), ta1_(v2), ta2_(I1620) { }
};

class Task874 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task874(std::shared_ptr<TATensor<double,4>> I1620, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1620), ta1_(Gamma503), ta2_(t2) { }
};

class Task875 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task875(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1623)
   : ta0_(I194), ta1_(v2), ta2_(I1623) { }
};

class Task876 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task876(std::shared_ptr<TATensor<double,4>> I1623, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1623), ta1_(Gamma51), ta2_(t2) { }
};

class Task877 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task877(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> Gamma562, std::shared_ptr<TATensor<double,4>> I1705)
   : ta0_(I194), ta1_(Gamma562), ta2_(I1705) { }
};

class Task878 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task878(std::shared_ptr<TATensor<double,4>> I1705, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1705), ta1_(t2) { }
};

class Task879 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task879(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> Gamma564, std::shared_ptr<TATensor<double,4>> I1709)
   : ta0_(I194), ta1_(Gamma564), ta2_(I1709) { }
};

class Task880 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task880(std::shared_ptr<TATensor<double,4>> I1709, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1709), ta1_(t2) { }
};

class Task881 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task881(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I239)
   : ta0_(proj), ta1_(I239) { }
};

class Task882 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task882(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I240)
   : ta0_(I239), ta1_(h1), ta2_(I240) { }
};

class Task883 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task883(std::shared_ptr<TATensor<double,4>> I240, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I240), ta1_(Gamma50), ta2_(t2) { }
};

class Task884 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task884(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I243)
   : ta0_(I239), ta1_(Gamma51), ta2_(I243) { }
};

class Task885 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task885(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I243), ta1_(t2), ta2_(h1) { }
};

class Task886 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task886(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I243), ta1_(t2), ta2_(h1) { }
};

class Task887 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task887(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task888 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task888(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task889 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task889(std::shared_ptr<TATensor<double,4>> I243, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I243), ta1_(t2), ta2_(v2) { }
};

class Task890 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task890(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1626)
   : ta0_(I239), ta1_(t2), ta2_(I1626) { }
};

class Task891 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task891(std::shared_ptr<TATensor<double,6>> I1626, std::shared_ptr<TATensor<double,6>> Gamma531, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1626), ta1_(Gamma531), ta2_(v2) { }
};

class Task892 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task892(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1629)
   : ta0_(I239), ta1_(t2), ta2_(I1629) { }
};

class Task893 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task893(std::shared_ptr<TATensor<double,6>> I1629, std::shared_ptr<TATensor<double,6>> Gamma532, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1629), ta1_(Gamma532), ta2_(v2) { }
};

class Task894 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task894(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I1632)
   : ta0_(I239), ta1_(t2), ta2_(I1632) { }
};

class Task895 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task895(std::shared_ptr<TATensor<double,6>> I1632, std::shared_ptr<TATensor<double,8>> Gamma533, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1632), ta1_(Gamma533), ta2_(v2) { }
};

class Task896 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task896(std::shared_ptr<TATensor<double,6>> I1632, std::shared_ptr<TATensor<double,8>> Gamma349, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1632), ta1_(Gamma349), ta2_(v2) { }
};

class Task897 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task897(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1638)
   : ta0_(I239), ta1_(v2), ta2_(I1638) { }
};

class Task898 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task898(std::shared_ptr<TATensor<double,4>> I1638, std::shared_ptr<TATensor<double,6>> Gamma471, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1638), ta1_(Gamma471), ta2_(t2) { }
};

class Task899 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task899(std::shared_ptr<TATensor<double,4>> I239, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1644)
   : ta0_(I239), ta1_(t2), ta2_(I1644) { }
};


}
}
}
#endif
#endif

