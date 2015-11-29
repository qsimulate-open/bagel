//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks3.h
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

#ifndef __SRC_SMITH_MRCI_TASKS3_H
#define __SRC_SMITH_MRCI_TASKS3_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task100 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<8,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,8>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<8,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task100(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task100() {}
};

class Task101 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task101(std::array<std::shared_ptr<Tensor>,2> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task101() {}
};

class Task102 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,5>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task102(std::array<std::shared_ptr<Tensor>,6> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task102() {}
};

class Task103 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<8,6> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,8>& block, const std::array<std::shared_ptr<const Tensor>,6>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<8,6>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task103(std::array<std::shared_ptr<Tensor>,7> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task103() {}
};

class Task104 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task104(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task104() {}
};

class Task105 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,3> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,3>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,3>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task105(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task105() {}
};

class Task106 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,3> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,3>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,3>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task106(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task106() {}
};

class Task107 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<8,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,8>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<8,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task107(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task107() {}
};

class Task108 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task108(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task109 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task109(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I0)
   : ta0_(proj), ta1_(I0) { }
};

class Task110 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task110(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> I1)
   : ta0_(I0), ta1_(Gamma0), ta2_(I1) { }
};

class Task111 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task111(std::shared_ptr<TATensor<double,4>> I1, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I1), ta1_(t2), ta2_(h1) { }
};

class Task112 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task112(std::shared_ptr<TATensor<double,4>> I1, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I280)
   : ta0_(I1), ta1_(t2), ta2_(I280) { }
};

class Task113 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task113(std::shared_ptr<TATensor<double,4>> I280, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I280), ta1_(v2) { }
};

class Task114 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task114(std::shared_ptr<TATensor<double,4>> I1, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1), ta1_(t2), ta2_(v2) { }
};

class Task115 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task115(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I4)
   : ta0_(I0), ta1_(h1), ta2_(I4) { }
};

class Task116 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task116(std::shared_ptr<TATensor<double,4>> I4, std::shared_ptr<TATensor<double,6>> Gamma1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I4), ta1_(Gamma1), ta2_(t2) { }
};

class Task117 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task117(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> Gamma2, std::shared_ptr<TATensor<double,4>> I7)
   : ta0_(I0), ta1_(Gamma2), ta2_(I7) { }
};

class Task118 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task118(std::shared_ptr<TATensor<double,4>> I7, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I7), ta1_(t2), ta2_(h1) { }
};

class Task119 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task119(std::shared_ptr<TATensor<double,4>> I7, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I7), ta1_(t2), ta2_(v2) { }
};

class Task120 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task120(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,6>> Gamma80, std::shared_ptr<TATensor<double,6>> I249)
   : ta0_(I0), ta1_(Gamma80), ta2_(I249) { }
};

class Task121 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task121(std::shared_ptr<TATensor<double,6>> I249, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I250)
   : ta0_(I249), ta1_(t2), ta2_(I250) { }
};

class Task122 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task122(std::shared_ptr<TATensor<double,4>> I250, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I250), ta1_(v2) { }
};

class Task123 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task123(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,6>> Gamma81, std::shared_ptr<TATensor<double,6>> I252)
   : ta0_(I0), ta1_(Gamma81), ta2_(I252) { }
};

class Task124 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task124(std::shared_ptr<TATensor<double,6>> I252, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I252), ta1_(t2), ta2_(v2) { }
};

class Task125 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task125(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,6>> Gamma82, std::shared_ptr<TATensor<double,6>> I255)
   : ta0_(I0), ta1_(Gamma82), ta2_(I255) { }
};

class Task126 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task126(std::shared_ptr<TATensor<double,6>> I255, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I255), ta1_(t2), ta2_(v2) { }
};

class Task127 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task127(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,6>> I264)
   : ta0_(I0), ta1_(t2), ta2_(I264) { }
};

class Task128 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task128(std::shared_ptr<TATensor<double,6>> I264, std::shared_ptr<TATensor<double,8>> Gamma85, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I264), ta1_(Gamma85), ta2_(v2) { }
};

class Task129 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task129(std::shared_ptr<TATensor<double,6>> I264, std::shared_ptr<TATensor<double,8>> Gamma86, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I264), ta1_(Gamma86), ta2_(v2) { }
};

class Task130 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task130(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I270)
   : ta0_(I0), ta1_(v2), ta2_(I270) { }
};

class Task131 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task131(std::shared_ptr<TATensor<double,4>> I270, std::shared_ptr<TATensor<double,6>> Gamma87, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I270), ta1_(Gamma87), ta2_(t2) { }
};

class Task132 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task132(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I273)
   : ta0_(I0), ta1_(t2), ta2_(I273) { }
};

class Task133 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task133(std::shared_ptr<TATensor<double,4>> I273, std::shared_ptr<TATensor<double,6>> Gamma88, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I273), ta1_(Gamma88), ta2_(v2) { }
};

class Task134 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task134(std::shared_ptr<TATensor<double,4>> I273, std::shared_ptr<TATensor<double,6>> Gamma89, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I273), ta1_(Gamma89), ta2_(v2) { }
};

class Task135 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task135(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,6>> Gamma94, std::shared_ptr<TATensor<double,6>> I291)
   : ta0_(I0), ta1_(Gamma94), ta2_(I291) { }
};

class Task136 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task136(std::shared_ptr<TATensor<double,6>> I291, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I291), ta1_(t2), ta2_(v2) { }
};

class Task137 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task137(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,6>> Gamma87, std::shared_ptr<TATensor<double,6>> I294)
   : ta0_(I0), ta1_(Gamma87), ta2_(I294) { }
};

class Task138 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task138(std::shared_ptr<TATensor<double,6>> I294, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I294), ta1_(t2), ta2_(v2) { }
};

class Task139 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task139(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I9)
   : ta0_(proj), ta1_(I9) { }
};

class Task140 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task140(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,6>> Gamma3, std::shared_ptr<TATensor<double,4>> I10)
   : ta0_(I9), ta1_(Gamma3), ta2_(I10) { }
};

class Task141 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task141(std::shared_ptr<TATensor<double,4>> I10, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I10), ta1_(t2), ta2_(h1) { }
};

class Task142 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task142(std::shared_ptr<TATensor<double,4>> I10, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I10), ta1_(t2), ta2_(v2) { }
};

class Task143 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task143(std::shared_ptr<TATensor<double,4>> I10, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I10), ta1_(t2), ta2_(v2) { }
};

class Task144 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task144(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I13)
   : ta0_(I9), ta1_(h1), ta2_(I13) { }
};

class Task145 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task145(std::shared_ptr<TATensor<double,4>> I13, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I13), ta1_(Gamma4), ta2_(t2) { }
};

class Task146 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task146(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,2>> I16)
   : ta0_(I9), ta1_(Gamma5), ta2_(I16) { }
};

class Task147 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task147(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I16), ta1_(t2), ta2_(h1) { }
};

class Task148 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task148(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I16), ta1_(t2), ta2_(h1) { }
};

class Task149 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task149(std::shared_ptr<TATensor<double,2>> I16, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I16), ta1_(t2), ta2_(v2) { }
};


}
}
}
#endif
#endif

