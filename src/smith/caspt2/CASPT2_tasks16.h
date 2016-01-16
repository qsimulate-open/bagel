//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks16.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS16_H
#define __SRC_SMITH_CASPT2_TASKS16_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task750 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task750(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma307, std::shared_ptr<TATensor<double,6>> I926)
   : ta0_(I698), ta1_(Gamma307), ta2_(I926) { }
};

class Task751 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task751(std::shared_ptr<TATensor<double,6>> I926, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I927)
   : ta0_(I926), ta1_(t2), ta2_(I927) { }
};

class Task752 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task752(std::shared_ptr<TATensor<double,4>> I927, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I927), ta1_(f1), ta2_(t2) { }
};

class Task753 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task753(std::shared_ptr<TATensor<double,6>> I926, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I939)
   : ta0_(I926), ta1_(t2), ta2_(I939) { }
};

class Task754 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task754(std::shared_ptr<TATensor<double,4>> I939, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I939), ta1_(t2), ta2_(f1) { }
};

class Task755 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task755(std::shared_ptr<TATensor<double,6>> I926, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1047)
   : ta0_(I926), ta1_(t2), ta2_(I1047) { }
};

class Task756 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task756(std::shared_ptr<TATensor<double,4>> I1047, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I1047), ta1_(t2), e0_(e) { }
};

class Task757 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task757(std::shared_ptr<TATensor<double,4>> I1047, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1047), ta1_(f1), ta2_(t2) { }
};

class Task758 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task758(std::shared_ptr<TATensor<double,6>> I926, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I926), ta1_(v2), ta2_(t2) { }
};

class Task759 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task759(std::shared_ptr<TATensor<double,6>> I926, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I926), ta1_(v2), ta2_(t2) { }
};

class Task760 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task760(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma308, std::shared_ptr<TATensor<double,4>> I930)
   : ta0_(I698), ta1_(Gamma308), ta2_(I930) { }
};

class Task761 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task761(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I931)
   : ta0_(I930), ta1_(t2), ta2_(I931) { }
};

class Task762 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task762(std::shared_ptr<TATensor<double,2>> I931, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I931), ta1_(t2), ta2_(f1) { }
};

class Task763 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task763(std::shared_ptr<TATensor<double,2>> I931, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I931), ta1_(t2), ta2_(f1) { }
};

class Task764 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task764(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I997)
   : ta0_(I930), ta1_(t2), ta2_(I997) { }
};

class Task765 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task765(std::shared_ptr<TATensor<double,2>> I997, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I997), ta1_(f1), ta2_(t2) { }
};

class Task766 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task766(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1001)
   : ta0_(I930), ta1_(t2), ta2_(I1001) { }
};

class Task767 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task767(std::shared_ptr<TATensor<double,2>> I1001, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1001), ta1_(f1), ta2_(t2) { }
};

class Task768 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task768(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1043)
   : ta0_(I930), ta1_(t2), ta2_(I1043) { }
};

class Task769 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task769(std::shared_ptr<TATensor<double,4>> I1043, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1043), ta1_(f1), ta2_(t2) { }
};

class Task770 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task770(std::shared_ptr<TATensor<double,4>> I1043, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1043), ta1_(f1), ta2_(t2) { }
};

class Task771 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task771(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1051)
   : ta0_(I930), ta1_(t2), ta2_(I1051) { }
};

class Task772 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task772(std::shared_ptr<TATensor<double,4>> I1051, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I1051), ta1_(t2), ta2_(f1) { }
};

class Task773 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task773(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I930), ta1_(t2), e0_(e) { }
};

class Task774 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task774(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I930), ta1_(v2), ta2_(t2) { }
};

class Task775 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task775(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I930), ta1_(v2), ta2_(t2) { }
};

class Task776 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task776(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I930), ta1_(h1), ta2_(t2) { }
};

class Task777 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task777(std::shared_ptr<TATensor<double,4>> I930, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I930), ta1_(h1), ta2_(t2) { }
};

class Task778 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,1>> ta1_;
    std::shared_ptr<TATensor<double,0>> ta2_;
    void compute_() override;
  public:
    Task778(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,1>> Gamma317, std::shared_ptr<TATensor<double,0>> I966)
   : ta0_(I698), ta1_(Gamma317), ta2_(I966) { }
};

class Task779 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task779(std::shared_ptr<TATensor<double,0>> I966, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I966), ta1_(t2) { }
};

class Task780 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task780(std::shared_ptr<TATensor<double,0>> I966, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I966), ta1_(t2) { }
};

class Task781 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task781(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,3>> Gamma329, std::shared_ptr<TATensor<double,2>> I1012)
   : ta0_(I698), ta1_(Gamma329), ta2_(I1012) { }
};

class Task782 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task782(std::shared_ptr<TATensor<double,2>> I1012, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1012), ta1_(t2) { }
};

class Task783 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task783(std::shared_ptr<TATensor<double,2>> I1012, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1012), ta1_(t2) { }
};

class Task784 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task784(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,5>> Gamma340, std::shared_ptr<TATensor<double,4>> I1054)
   : ta0_(I698), ta1_(Gamma340), ta2_(I1054) { }
};

class Task785 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task785(std::shared_ptr<TATensor<double,4>> I1054, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1054), ta1_(t2) { }
};

class Task786 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task786(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma355, std::shared_ptr<TATensor<double,6>> I1100)
   : ta0_(I698), ta1_(Gamma355), ta2_(I1100) { }
};

class Task787 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task787(std::shared_ptr<TATensor<double,6>> I1100, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1100), ta1_(v2), ta2_(t2) { }
};

class Task788 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task788(std::shared_ptr<TATensor<double,1>> I698, std::shared_ptr<TATensor<double,7>> Gamma373, std::shared_ptr<TATensor<double,6>> I1154)
   : ta0_(I698), ta1_(Gamma373), ta2_(I1154) { }
};

class Task789 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task789(std::shared_ptr<TATensor<double,6>> I1154, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1154), ta1_(v2), ta2_(t2) { }
};


}
}
}
#endif
#endif

