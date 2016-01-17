//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks16.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS16_H
#define __SRC_SMITH_RelMRCI_TASKS16_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task750 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task750(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1337)
   : ta0_(proj), ta1_(I1337) { }
};

class Task751 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task751(std::shared_ptr<TATensor<std::complex<double>,4>> I1337, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1337), ta1_(Gamma4), ta2_(t2) { }
};

class Task752 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task752(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1339)
   : ta0_(proj), ta1_(I1339) { }
};

class Task753 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task753(std::shared_ptr<TATensor<std::complex<double>,4>> I1339, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I1340)
   : ta0_(I1339), ta1_(Gamma11), ta2_(I1340) { }
};

class Task754 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task754(std::shared_ptr<TATensor<std::complex<double>,4>> I1340, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1340), ta1_(t2) { }
};

class Task755 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task755(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1343)
   : ta0_(proj), ta1_(I1343) { }
};

class Task756 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task756(std::shared_ptr<TATensor<std::complex<double>,4>> I1343, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1343), ta1_(Gamma24), ta2_(t2) { }
};

class Task757 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task757(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1345)
   : ta0_(proj), ta1_(I1345) { }
};

class Task758 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task758(std::shared_ptr<TATensor<std::complex<double>,4>> I1345, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1345), ta1_(Gamma32), ta2_(t2) { }
};

class Task759 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task759(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1347)
   : ta0_(proj), ta1_(I1347) { }
};

class Task760 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task760(std::shared_ptr<TATensor<std::complex<double>,4>> I1347, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1347), ta1_(t2) { }
};

class Task761 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task761(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1349)
   : ta0_(proj), ta1_(I1349) { }
};

class Task762 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task762(std::shared_ptr<TATensor<std::complex<double>,4>> I1349, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I1350)
   : ta0_(I1349), ta1_(Gamma27), ta2_(I1350) { }
};

class Task763 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task763(std::shared_ptr<TATensor<std::complex<double>,4>> I1350, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1350), ta1_(t2) { }
};

class Task764 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task764(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1353)
   : ta0_(proj), ta1_(I1353) { }
};

class Task765 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task765(std::shared_ptr<TATensor<std::complex<double>,4>> I1353, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1353), ta1_(Gamma33), ta2_(t2) { }
};


}
}
}
#endif
#endif

