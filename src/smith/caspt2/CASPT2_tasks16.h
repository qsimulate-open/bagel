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
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task750(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I1279)
   : ta0_(I915), ta1_(h1), ta2_(I1279) { }
};

class Task751 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task751(std::shared_ptr<TATensor<double,4>> I1279, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1279), ta1_(t2) { }
};

class Task752 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task752(std::shared_ptr<TATensor<double,2>> I915, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I1282)
   : ta0_(I915), ta1_(h1), ta2_(I1282) { }
};

class Task753 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task753(std::shared_ptr<TATensor<double,4>> I1282, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1282), ta1_(t2) { }
};

class Task754 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task754(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma321, std::shared_ptr<TATensor<double,6>> I965)
   : ta0_(I768), ta1_(Gamma321), ta2_(I965) { }
};

class Task755 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task755(std::shared_ptr<TATensor<double,6>> I965, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I966)
   : ta0_(I965), ta1_(t2), ta2_(I966) { }
};

class Task756 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task756(std::shared_ptr<TATensor<double,4>> I966, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I966), ta1_(f1), ta2_(t2) { }
};

class Task757 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task757(std::shared_ptr<TATensor<double,6>> I965, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I965), ta1_(v2), ta2_(t2) { }
};

class Task758 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task758(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma326, std::shared_ptr<TATensor<double,6>> I985)
   : ta0_(I768), ta1_(Gamma326), ta2_(I985) { }
};

class Task759 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task759(std::shared_ptr<TATensor<double,6>> I985, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I986)
   : ta0_(I985), ta1_(t2), ta2_(I986) { }
};

class Task760 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task760(std::shared_ptr<TATensor<double,4>> I986, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I986), ta1_(t2), ta2_(f1) { }
};

class Task761 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task761(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma327, std::shared_ptr<TATensor<double,6>> I989)
   : ta0_(I768), ta1_(Gamma327), ta2_(I989) { }
};

class Task762 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task762(std::shared_ptr<TATensor<double,6>> I989, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I990)
   : ta0_(I989), ta1_(t2), ta2_(I990) { }
};

class Task763 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task763(std::shared_ptr<TATensor<double,4>> I990, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I990), ta1_(t2), ta2_(f1) { }
};

class Task764 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task764(std::shared_ptr<TATensor<double,6>> I989, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I989), ta1_(v2), ta2_(t2) { }
};

class Task765 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task765(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma328, std::shared_ptr<TATensor<double,6>> I993)
   : ta0_(I768), ta1_(Gamma328), ta2_(I993) { }
};

class Task766 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task766(std::shared_ptr<TATensor<double,6>> I993, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I993), ta1_(t2) { }
};

class Task767 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,7>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task767(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,7>> Gamma329, std::shared_ptr<TATensor<double,6>> I996)
   : ta0_(I768), ta1_(Gamma329), ta2_(I996) { }
};

class Task768 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task768(std::shared_ptr<TATensor<double,6>> I996, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I997)
   : ta0_(I996), ta1_(t2), ta2_(I997) { }
};

class Task769 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task769(std::shared_ptr<TATensor<double,4>> I997, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I997), ta1_(f1), ta2_(t2) { }
};

class Task770 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task770(std::shared_ptr<TATensor<double,6>> I996, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1009)
   : ta0_(I996), ta1_(t2), ta2_(I1009) { }
};

class Task771 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task771(std::shared_ptr<TATensor<double,4>> I1009, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I1009), ta1_(t2), ta2_(f1) { }
};

class Task772 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task772(std::shared_ptr<TATensor<double,6>> I996, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1117)
   : ta0_(I996), ta1_(t2), ta2_(I1117) { }
};

class Task773 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task773(std::shared_ptr<TATensor<double,4>> I1117, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I1117), ta1_(t2), e0_(e) { }
};

class Task774 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task774(std::shared_ptr<TATensor<double,4>> I1117, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1117), ta1_(f1), ta2_(t2) { }
};

class Task775 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task775(std::shared_ptr<TATensor<double,6>> I996, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I996), ta1_(v2), ta2_(t2) { }
};

class Task776 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task776(std::shared_ptr<TATensor<double,6>> I996, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I996), ta1_(v2), ta2_(t2) { }
};

class Task777 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task777(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma330, std::shared_ptr<TATensor<double,4>> I1000)
   : ta0_(I768), ta1_(Gamma330), ta2_(I1000) { }
};

class Task778 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task778(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1001)
   : ta0_(I1000), ta1_(t2), ta2_(I1001) { }
};

class Task779 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task779(std::shared_ptr<TATensor<double,2>> I1001, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I1001), ta1_(t2), ta2_(f1) { }
};

class Task780 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task780(std::shared_ptr<TATensor<double,2>> I1001, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I1001), ta1_(t2), ta2_(f1) { }
};

class Task781 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task781(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1067)
   : ta0_(I1000), ta1_(t2), ta2_(I1067) { }
};

class Task782 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task782(std::shared_ptr<TATensor<double,2>> I1067, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1067), ta1_(f1), ta2_(t2) { }
};

class Task783 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task783(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I1071)
   : ta0_(I1000), ta1_(t2), ta2_(I1071) { }
};

class Task784 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task784(std::shared_ptr<TATensor<double,2>> I1071, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1071), ta1_(f1), ta2_(t2) { }
};

class Task785 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task785(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1113)
   : ta0_(I1000), ta1_(t2), ta2_(I1113) { }
};

class Task786 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task786(std::shared_ptr<TATensor<double,4>> I1113, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1113), ta1_(f1), ta2_(t2) { }
};

class Task787 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task787(std::shared_ptr<TATensor<double,4>> I1113, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1113), ta1_(f1), ta2_(t2) { }
};

class Task788 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task788(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1121)
   : ta0_(I1000), ta1_(t2), ta2_(I1121) { }
};

class Task789 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task789(std::shared_ptr<TATensor<double,4>> I1121, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I1121), ta1_(t2), ta2_(f1) { }
};

class Task790 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task790(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I1000), ta1_(t2), e0_(e) { }
};

class Task791 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task791(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1000), ta1_(v2), ta2_(t2) { }
};

class Task792 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task792(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1000), ta1_(v2), ta2_(t2) { }
};

class Task793 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task793(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1000), ta1_(h1), ta2_(t2) { }
};

class Task794 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task794(std::shared_ptr<TATensor<double,4>> I1000, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1000), ta1_(h1), ta2_(t2) { }
};

class Task795 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,1>> ta1_;
    std::shared_ptr<TATensor<double,0>> ta2_;
    void compute_() override;
  public:
    Task795(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,1>> Gamma339, std::shared_ptr<TATensor<double,0>> I1036)
   : ta0_(I768), ta1_(Gamma339), ta2_(I1036) { }
};

class Task796 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task796(std::shared_ptr<TATensor<double,0>> I1036, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1036), ta1_(t2) { }
};

class Task797 : public Task {
  protected:
    std::shared_ptr<TATensor<double,0>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task797(std::shared_ptr<TATensor<double,0>> I1036, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1036), ta1_(t2) { }
};

class Task798 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,3>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task798(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,3>> Gamma351, std::shared_ptr<TATensor<double,2>> I1082)
   : ta0_(I768), ta1_(Gamma351), ta2_(I1082) { }
};

class Task799 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task799(std::shared_ptr<TATensor<double,2>> I1082, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1082), ta1_(t2) { }
};


}
}
}
#endif
#endif

