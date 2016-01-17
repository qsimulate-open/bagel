//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks17.h
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

#ifndef __SRC_SMITH_MRCI_TASKS17_H
#define __SRC_SMITH_MRCI_TASKS17_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task800 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task800(std::shared_ptr<TATensor<double,4>> I1394, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1394), ta1_(Gamma27), ta2_(t2) { }
};

class Task801 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task801(std::shared_ptr<TATensor<double,4>> I1394, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1394), ta1_(Gamma29), ta2_(t2) { }
};

class Task802 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task802(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1397)
   : ta0_(I194), ta1_(v2), ta2_(I1397) { }
};

class Task803 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task803(std::shared_ptr<TATensor<double,4>> I1397, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1398)
   : ta0_(I1397), ta1_(Gamma29), ta2_(I1398) { }
};

class Task804 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task804(std::shared_ptr<TATensor<double,4>> I1398, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1398), ta1_(t2) { }
};

class Task805 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task805(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1400)
   : ta0_(I194), ta1_(t2), ta2_(I1400) { }
};

class Task806 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task806(std::shared_ptr<TATensor<double,4>> I1400, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1400), ta1_(Gamma49), ta2_(v2) { }
};

class Task807 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task807(std::shared_ptr<TATensor<double,4>> I1400, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1400), ta1_(Gamma240), ta2_(v2) { }
};

class Task808 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task808(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1403)
   : ta0_(I194), ta1_(t2), ta2_(I1403) { }
};

class Task809 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task809(std::shared_ptr<TATensor<double,4>> I1403, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1403), ta1_(Gamma49), ta2_(v2) { }
};

class Task810 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task810(std::shared_ptr<TATensor<double,4>> I1403, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1403), ta1_(Gamma240), ta2_(v2) { }
};

class Task811 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task811(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1430)
   : ta0_(I194), ta1_(v2), ta2_(I1430) { }
};

class Task812 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task812(std::shared_ptr<TATensor<double,4>> I1430, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1430), ta1_(Gamma50), ta2_(t2) { }
};

class Task813 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task813(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1433)
   : ta0_(I194), ta1_(v2), ta2_(I1433) { }
};

class Task814 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task814(std::shared_ptr<TATensor<double,4>> I1433, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1433), ta1_(Gamma50), ta2_(t2) { }
};

class Task815 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task815(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1436)
   : ta0_(I194), ta1_(v2), ta2_(I1436) { }
};

class Task816 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task816(std::shared_ptr<TATensor<double,4>> I1436, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1436), ta1_(Gamma252), ta2_(t2) { }
};

class Task817 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task817(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1439)
   : ta0_(I194), ta1_(v2), ta2_(I1439) { }
};

class Task818 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task818(std::shared_ptr<TATensor<double,4>> I1439, std::shared_ptr<TATensor<double,6>> Gamma31, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1439), ta1_(Gamma31), ta2_(t2) { }
};

class Task819 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task819(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1442)
   : ta0_(I194), ta1_(v2), ta2_(I1442) { }
};

class Task820 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task820(std::shared_ptr<TATensor<double,4>> I1442, std::shared_ptr<TATensor<double,6>> Gamma471, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1442), ta1_(Gamma471), ta2_(t2) { }
};

class Task821 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task821(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1445)
   : ta0_(I194), ta1_(v2), ta2_(I1445) { }
};

class Task822 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task822(std::shared_ptr<TATensor<double,4>> I1445, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1445), ta1_(Gamma50), ta2_(t2) { }
};

class Task823 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task823(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1448)
   : ta0_(I194), ta1_(v2), ta2_(I1448) { }
};

class Task824 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task824(std::shared_ptr<TATensor<double,4>> I1448, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1448), ta1_(Gamma50), ta2_(t2) { }
};

class Task825 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task825(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1451)
   : ta0_(I194), ta1_(v2), ta2_(I1451) { }
};

class Task826 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task826(std::shared_ptr<TATensor<double,4>> I1451, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1451), ta1_(Gamma50), ta2_(t2) { }
};

class Task827 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task827(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1454)
   : ta0_(I194), ta1_(v2), ta2_(I1454) { }
};

class Task828 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task828(std::shared_ptr<TATensor<double,2>> I1454, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1454), ta1_(Gamma51), ta2_(t2) { }
};

class Task829 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task829(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1457)
   : ta0_(I194), ta1_(v2), ta2_(I1457) { }
};

class Task830 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task830(std::shared_ptr<TATensor<double,2>> I1457, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1457), ta1_(Gamma51), ta2_(t2) { }
};

class Task831 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task831(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1472)
   : ta0_(I194), ta1_(t2), ta2_(I1472) { }
};

class Task832 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task832(std::shared_ptr<TATensor<double,4>> I1472, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1472), ta1_(Gamma32), ta2_(v2) { }
};

class Task833 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task833(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1475)
   : ta0_(I194), ta1_(t2), ta2_(I1475) { }
};

class Task834 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task834(std::shared_ptr<TATensor<double,4>> I1475, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1475), ta1_(Gamma32), ta2_(v2) { }
};

class Task835 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task835(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1502)
   : ta0_(I194), ta1_(t2), ta2_(I1502) { }
};

class Task836 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task836(std::shared_ptr<TATensor<double,4>> I1502, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1503)
   : ta0_(I1502), ta1_(Gamma51), ta2_(I1503) { }
};

class Task837 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task837(std::shared_ptr<TATensor<double,4>> I1503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1503), ta1_(v2) { }
};

class Task838 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task838(std::shared_ptr<TATensor<double,4>> I1502, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1502), ta1_(Gamma29), ta2_(v2) { }
};

class Task839 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task839(std::shared_ptr<TATensor<double,4>> I1502, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1502), ta1_(Gamma503), ta2_(v2) { }
};

class Task840 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task840(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1505)
   : ta0_(I194), ta1_(t2), ta2_(I1505) { }
};

class Task841 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task841(std::shared_ptr<TATensor<double,4>> I1505, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1506)
   : ta0_(I1505), ta1_(Gamma51), ta2_(I1506) { }
};

class Task842 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task842(std::shared_ptr<TATensor<double,4>> I1506, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1506), ta1_(v2) { }
};

class Task843 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task843(std::shared_ptr<TATensor<double,4>> I1505, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1505), ta1_(Gamma27), ta2_(v2) { }
};

class Task844 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task844(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1508)
   : ta0_(I194), ta1_(t2), ta2_(I1508) { }
};

class Task845 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task845(std::shared_ptr<TATensor<double,4>> I1508, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1509)
   : ta0_(I1508), ta1_(Gamma51), ta2_(I1509) { }
};

class Task846 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task846(std::shared_ptr<TATensor<double,4>> I1509, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1509), ta1_(v2) { }
};

class Task847 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task847(std::shared_ptr<TATensor<double,4>> I1508, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1508), ta1_(Gamma29), ta2_(v2) { }
};

class Task848 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task848(std::shared_ptr<TATensor<double,4>> I1508, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1508), ta1_(Gamma503), ta2_(v2) { }
};

class Task849 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task849(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1511)
   : ta0_(I194), ta1_(t2), ta2_(I1511) { }
};


}
}
}
#endif
#endif

