//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks16.h
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

#ifndef __SRC_SMITH_MRCI_TASKS16_H
#define __SRC_SMITH_MRCI_TASKS16_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task750 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task750(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I219), ta1_(t2), ta2_(h1) { }
};

class Task751 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task751(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I219), ta1_(t2), ta2_(h1) { }
};

class Task752 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task752(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I219), ta1_(t2), ta2_(h1) { }
};

class Task753 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task753(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task754 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task754(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task755 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task755(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task756 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task756(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task757 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task757(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task758 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task758(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task759 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task759(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task760 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task760(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task761 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task761(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task762 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task762(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task763 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task763(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task764 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task764(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task765 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task765(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task766 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task766(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task767 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task767(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task768 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task768(std::shared_ptr<TATensor<double,4>> I219, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I219), ta1_(t2), ta2_(v2) { }
};

class Task769 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task769(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I237)
   : ta0_(I194), ta1_(h1), ta2_(I237) { }
};

class Task770 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task770(std::shared_ptr<TATensor<double,4>> I237, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I237), ta1_(Gamma51), ta2_(t2) { }
};

class Task771 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task771(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1355)
   : ta0_(I194), ta1_(v2), ta2_(I1355) { }
};

class Task772 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task772(std::shared_ptr<TATensor<double,4>> I1355, std::shared_ptr<TATensor<double,6>> Gamma24, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1355), ta1_(Gamma24), ta2_(t2) { }
};

class Task773 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task773(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1358)
   : ta0_(I194), ta1_(t2), ta2_(I1358) { }
};

class Task774 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task774(std::shared_ptr<TATensor<double,4>> I1358, std::shared_ptr<TATensor<double,4>> Gamma25, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1358), ta1_(Gamma25), ta2_(v2) { }
};

class Task775 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task775(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1361)
   : ta0_(I194), ta1_(t2), ta2_(I1361) { }
};

class Task776 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task776(std::shared_ptr<TATensor<double,4>> I1361, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1361), ta1_(Gamma5), ta2_(v2) { }
};

class Task777 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task777(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1364)
   : ta0_(I194), ta1_(t2), ta2_(I1364) { }
};

class Task778 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task778(std::shared_ptr<TATensor<double,4>> I1364, std::shared_ptr<TATensor<double,4>> Gamma25, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1364), ta1_(Gamma25), ta2_(v2) { }
};

class Task779 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task779(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1367)
   : ta0_(I194), ta1_(t2), ta2_(I1367) { }
};

class Task780 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task780(std::shared_ptr<TATensor<double,4>> I1367, std::shared_ptr<TATensor<double,4>> Gamma25, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1367), ta1_(Gamma25), ta2_(v2) { }
};

class Task781 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task781(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1370)
   : ta0_(I194), ta1_(t2), ta2_(I1370) { }
};

class Task782 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task782(std::shared_ptr<TATensor<double,4>> I1370, std::shared_ptr<TATensor<double,6>> Gamma49, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1370), ta1_(Gamma49), ta2_(v2) { }
};

class Task783 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task783(std::shared_ptr<TATensor<double,4>> I1370, std::shared_ptr<TATensor<double,6>> Gamma240, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1370), ta1_(Gamma240), ta2_(v2) { }
};

class Task784 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task784(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1373)
   : ta0_(I194), ta1_(t2), ta2_(I1373) { }
};

class Task785 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task785(std::shared_ptr<TATensor<double,4>> I1373, std::shared_ptr<TATensor<double,6>> Gamma48, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1373), ta1_(Gamma48), ta2_(v2) { }
};

class Task786 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task786(std::shared_ptr<TATensor<double,4>> I1373, std::shared_ptr<TATensor<double,6>> Gamma230, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1373), ta1_(Gamma230), ta2_(v2) { }
};

class Task787 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task787(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1382)
   : ta0_(I194), ta1_(v2), ta2_(I1382) { }
};

class Task788 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task788(std::shared_ptr<TATensor<double,4>> I1382, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1382), ta1_(Gamma27), ta2_(t2) { }
};

class Task789 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task789(std::shared_ptr<TATensor<double,4>> I1382, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1382), ta1_(Gamma29), ta2_(t2) { }
};

class Task790 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task790(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1385)
   : ta0_(I194), ta1_(v2), ta2_(I1385) { }
};

class Task791 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task791(std::shared_ptr<TATensor<double,4>> I1385, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1385), ta1_(Gamma27), ta2_(t2) { }
};

class Task792 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task792(std::shared_ptr<TATensor<double,4>> I1385, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1385), ta1_(Gamma29), ta2_(t2) { }
};

class Task793 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task793(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1388)
   : ta0_(I194), ta1_(v2), ta2_(I1388) { }
};

class Task794 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task794(std::shared_ptr<TATensor<double,4>> I1388, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1389)
   : ta0_(I1388), ta1_(Gamma29), ta2_(I1389) { }
};

class Task795 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task795(std::shared_ptr<TATensor<double,4>> I1389, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1389), ta1_(t2) { }
};

class Task796 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task796(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1391)
   : ta0_(I194), ta1_(v2), ta2_(I1391) { }
};

class Task797 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task797(std::shared_ptr<TATensor<double,4>> I1391, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1391), ta1_(Gamma27), ta2_(t2) { }
};

class Task798 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task798(std::shared_ptr<TATensor<double,4>> I1391, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1391), ta1_(Gamma29), ta2_(t2) { }
};

class Task799 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task799(std::shared_ptr<TATensor<double,4>> I194, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1394)
   : ta0_(I194), ta1_(v2), ta2_(I1394) { }
};


}
}
}
#endif
#endif

