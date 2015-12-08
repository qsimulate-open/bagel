//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks20.h
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

#ifndef __SRC_SMITH_MRCI_TASKS20_H
#define __SRC_SMITH_MRCI_TASKS20_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task950 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task950(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1727)
   : ta0_(proj), ta1_(I1727) { }
};

class Task951 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task951(std::shared_ptr<TATensor<double,4>> I1727, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1727), ta1_(Gamma0), ta2_(v2) { }
};

class Task952 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task952(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1733)
   : ta0_(proj), ta1_(I1733) { }
};

class Task953 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task953(std::shared_ptr<TATensor<double,4>> I1733, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I1734)
   : ta0_(I1733), ta1_(Gamma12), ta2_(I1734) { }
};

class Task954 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task954(std::shared_ptr<TATensor<double,4>> I1734, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1734), ta1_(v2) { }
};

class Task955 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task955(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1757)
   : ta0_(proj), ta1_(I1757) { }
};

class Task956 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task956(std::shared_ptr<TATensor<double,4>> I1757, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1757), ta1_(v2) { }
};

class Task957 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task957(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1759)
   : ta0_(proj), ta1_(I1759) { }
};

class Task958 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task958(std::shared_ptr<TATensor<double,4>> I1759, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1760)
   : ta0_(I1759), ta1_(Gamma32), ta2_(I1760) { }
};

class Task959 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task959(std::shared_ptr<TATensor<double,4>> I1760, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1760), ta1_(v2) { }
};

class Task960 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task960(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1763)
   : ta0_(proj), ta1_(I1763) { }
};

class Task961 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task961(std::shared_ptr<TATensor<double,4>> I1763, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1763), ta1_(Gamma51), ta2_(v2) { }
};

class Task962 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task962(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task963 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task963(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1765)
   : ta0_(proj), ta1_(I1765) { }
};

class Task964 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task964(std::shared_ptr<TATensor<double,4>> I1765, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1765), ta1_(Gamma0), ta2_(t2) { }
};

class Task965 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task965(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1767)
   : ta0_(proj), ta1_(I1767) { }
};

class Task966 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task966(std::shared_ptr<TATensor<double,4>> I1767, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1767), ta1_(Gamma4), ta2_(t2) { }
};

class Task967 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task967(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1769)
   : ta0_(proj), ta1_(I1769) { }
};

class Task968 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task968(std::shared_ptr<TATensor<double,4>> I1769, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I1770)
   : ta0_(I1769), ta1_(Gamma12), ta2_(I1770) { }
};

class Task969 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task969(std::shared_ptr<TATensor<double,4>> I1770, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1770), ta1_(t2) { }
};

class Task970 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task970(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1773)
   : ta0_(proj), ta1_(I1773) { }
};

class Task971 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task971(std::shared_ptr<TATensor<double,4>> I1773, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1773), ta1_(Gamma27), ta2_(t2) { }
};

class Task972 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task972(std::shared_ptr<TATensor<double,4>> I1773, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1773), ta1_(Gamma29), ta2_(t2) { }
};

class Task973 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task973(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1777)
   : ta0_(proj), ta1_(I1777) { }
};

class Task974 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task974(std::shared_ptr<TATensor<double,4>> I1777, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1778)
   : ta0_(I1777), ta1_(Gamma29), ta2_(I1778) { }
};

class Task975 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task975(std::shared_ptr<TATensor<double,4>> I1778, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1778), ta1_(t2) { }
};

class Task976 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task976(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1781)
   : ta0_(proj), ta1_(I1781) { }
};

class Task977 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task977(std::shared_ptr<TATensor<double,4>> I1781, std::shared_ptr<TATensor<double,6>> Gamma50, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1781), ta1_(Gamma50), ta2_(t2) { }
};

class Task978 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task978(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1783)
   : ta0_(proj), ta1_(I1783) { }
};

class Task979 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task979(std::shared_ptr<TATensor<double,4>> I1783, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1783), ta1_(t2) { }
};

class Task980 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task980(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1785)
   : ta0_(proj), ta1_(I1785) { }
};

class Task981 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task981(std::shared_ptr<TATensor<double,4>> I1785, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1786)
   : ta0_(I1785), ta1_(Gamma32), ta2_(I1786) { }
};

class Task982 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task982(std::shared_ptr<TATensor<double,4>> I1786, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1786), ta1_(t2) { }
};

class Task983 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task983(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I1789)
   : ta0_(proj), ta1_(I1789) { }
};

class Task984 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task984(std::shared_ptr<TATensor<double,4>> I1789, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1789), ta1_(Gamma51), ta2_(t2) { }
};


}
}
}
#endif
#endif

