//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks17.h
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
    Task800(std::shared_ptr<TATensor<double,4>> I1395, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1395), ta1_(Gamma27), ta2_(t2) { }
};

class Task801 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task801(std::shared_ptr<TATensor<double,4>> I1395, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1395), ta1_(Gamma29), ta2_(t2) { }
};

class Task802 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task802(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1398)
   : ta0_(I194), ta1_(v2), ta2_(I1398) { }
};

class Task803 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task803(std::shared_ptr<TATensor<double,4>> I1398, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1398), ta1_(Gamma27), ta2_(t2) { }
};

class Task804 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task804(std::shared_ptr<TATensor<double,4>> I1398, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1398), ta1_(Gamma29), ta2_(t2) { }
};

class Task805 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task805(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1401)
   : ta0_(I194), ta1_(v2), ta2_(I1401) { }
};

class Task806 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task806(std::shared_ptr<TATensor<double,4>> I1401, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1402)
   : ta0_(I1401), ta1_(Gamma29), ta2_(I1402) { }
};

class Task807 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task807(std::shared_ptr<TATensor<double,4>> I1402, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1402), ta1_(t2) { }
};

class Task808 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task808(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1404)
   : ta0_(I194), ta1_(t2), ta2_(I1404) { }
};

class Task809 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task809(std::shared_ptr<TATensor<double,4>> I1404, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1404), ta1_(Gamma49), ta2_(v2) { }
};

class Task810 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task810(std::shared_ptr<TATensor<double,4>> I1404, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1404), ta1_(Gamma240), ta2_(v2) { }
};

class Task811 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task811(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1407)
   : ta0_(I194), ta1_(t2), ta2_(I1407) { }
};

class Task812 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task812(std::shared_ptr<TATensor<double,4>> I1407, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1407), ta1_(Gamma49), ta2_(v2) { }
};

class Task813 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task813(std::shared_ptr<TATensor<double,4>> I1407, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1407), ta1_(Gamma240), ta2_(v2) { }
};

class Task814 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task814(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1434)
   : ta0_(I194), ta1_(v2), ta2_(I1434) { }
};

class Task815 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task815(std::shared_ptr<TATensor<double,4>> I1434, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1434), ta1_(Gamma50), ta2_(t2) { }
};

class Task816 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task816(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1437)
   : ta0_(I194), ta1_(v2), ta2_(I1437) { }
};

class Task817 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task817(std::shared_ptr<TATensor<double,4>> I1437, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1437), ta1_(Gamma50), ta2_(t2) { }
};

class Task818 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task818(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1440)
   : ta0_(I194), ta1_(v2), ta2_(I1440) { }
};

class Task819 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task819(std::shared_ptr<TATensor<double,4>> I1440, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1440), ta1_(Gamma252), ta2_(t2) { }
};

class Task820 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task820(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1443)
   : ta0_(I194), ta1_(v2), ta2_(I1443) { }
};

class Task821 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task821(std::shared_ptr<TATensor<double,4>> I1443, std::shared_ptr<TATensor<double,6>> Gamma31, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1443), ta1_(Gamma31), ta2_(t2) { }
};

class Task822 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task822(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1446)
   : ta0_(I194), ta1_(v2), ta2_(I1446) { }
};

class Task823 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task823(std::shared_ptr<TATensor<double,4>> I1446, std::shared_ptr<TATensor<double,6>> Gamma471, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1446), ta1_(Gamma471), ta2_(t2) { }
};

class Task824 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task824(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1449)
   : ta0_(I194), ta1_(v2), ta2_(I1449) { }
};

class Task825 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task825(std::shared_ptr<TATensor<double,4>> I1449, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1449), ta1_(Gamma50), ta2_(t2) { }
};

class Task826 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task826(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1452)
   : ta0_(I194), ta1_(v2), ta2_(I1452) { }
};

class Task827 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task827(std::shared_ptr<TATensor<double,4>> I1452, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1452), ta1_(Gamma50), ta2_(t2) { }
};

class Task828 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task828(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1455)
   : ta0_(I194), ta1_(v2), ta2_(I1455) { }
};

class Task829 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task829(std::shared_ptr<TATensor<double,4>> I1455, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1455), ta1_(Gamma50), ta2_(t2) { }
};

class Task830 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task830(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1458)
   : ta0_(I194), ta1_(v2), ta2_(I1458) { }
};

class Task831 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task831(std::shared_ptr<TATensor<double,2>> I1458, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1458), ta1_(Gamma51), ta2_(t2) { }
};

class Task832 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task832(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1461)
   : ta0_(I194), ta1_(v2), ta2_(I1461) { }
};

class Task833 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task833(std::shared_ptr<TATensor<double,2>> I1461, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1461), ta1_(Gamma51), ta2_(t2) { }
};

class Task834 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task834(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1476)
   : ta0_(I194), ta1_(t2), ta2_(I1476) { }
};

class Task835 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task835(std::shared_ptr<TATensor<double,4>> I1476, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1476), ta1_(Gamma32), ta2_(v2) { }
};

class Task836 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task836(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1479)
   : ta0_(I194), ta1_(t2), ta2_(I1479) { }
};

class Task837 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task837(std::shared_ptr<TATensor<double,4>> I1479, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1479), ta1_(Gamma32), ta2_(v2) { }
};

class Task838 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task838(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1506)
   : ta0_(I194), ta1_(t2), ta2_(I1506) { }
};

class Task839 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task839(std::shared_ptr<TATensor<double,4>> I1506, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1507)
   : ta0_(I1506), ta1_(Gamma51), ta2_(I1507) { }
};

class Task840 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task840(std::shared_ptr<TATensor<double,4>> I1507, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1507), ta1_(v2) { }
};

class Task841 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task841(std::shared_ptr<TATensor<double,4>> I1506, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1506), ta1_(Gamma29), ta2_(v2) { }
};

class Task842 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task842(std::shared_ptr<TATensor<double,4>> I1506, std::shared_ptr<TATensor<double,4>> Gamma503, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1506), ta1_(Gamma503), ta2_(v2) { }
};

class Task843 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task843(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1509)
   : ta0_(I194), ta1_(t2), ta2_(I1509) { }
};

class Task844 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task844(std::shared_ptr<TATensor<double,4>> I1509, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1510)
   : ta0_(I1509), ta1_(Gamma51), ta2_(I1510) { }
};

class Task845 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task845(std::shared_ptr<TATensor<double,4>> I1510, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1510), ta1_(v2) { }
};

class Task846 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task846(std::shared_ptr<TATensor<double,4>> I1509, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1509), ta1_(Gamma27), ta2_(v2) { }
};

class Task847 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task847(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1512)
   : ta0_(I194), ta1_(t2), ta2_(I1512) { }
};

class Task848 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task848(std::shared_ptr<TATensor<double,4>> I1512, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> I1513)
   : ta0_(I1512), ta1_(Gamma51), ta2_(I1513) { }
};

class Task849 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task849(std::shared_ptr<TATensor<double,4>> I1513, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1513), ta1_(v2) { }
};


}
}
}
#endif
#endif

