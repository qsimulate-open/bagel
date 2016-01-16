//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks2.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS2_H
#define __SRC_SMITH_CASPT2_TASKS2_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task50 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,3> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,3>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,3>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task50(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task50() {}
};

class Task51 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task51(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task51() {}
};

class Task52 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task52(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task52() {}
};

class Task53 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task53(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task53() {}
};

class Task54 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task54(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task54() {}
};

class Task55 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task55(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task55() {}
};

class Task56 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task56(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task56() {}
};

class Task57 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<3,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task57(std::array<std::shared_ptr<Tensor>,2> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task57() {}
};

class Task58 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task58(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task58() {}
};

class Task59 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task59(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task59() {}
};

class Task60 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task60(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task60() {}
};

class Task61 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<9,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,9>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<9,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task61(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task61() {}
};

class Task62 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task62(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task62() {}
};

class Task63 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task63(std::array<std::shared_ptr<Tensor>,2> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task63() {}
};

class Task64 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<3,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task64(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task64() {}
};

class Task65 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task65(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task65() {}
};

class Task66 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task66(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task66() {}
};

class Task67 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,3> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,3>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,3>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task67(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task67() {}
};

class Task68 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,3> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,3>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,3>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task68(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task68() {}
};

class Task69 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task69(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task70 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task70(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I0)
   : ta0_(proj), ta1_(I0) { }
};

class Task71 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task71(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> Gamma0, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I0), ta1_(Gamma0), ta2_(t2) { }
};

class Task72 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task72(std::shared_ptr<TATensor<double,4>> I0, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I0), ta1_(Gamma92), ta2_(t2), e0_(e) { }
};

class Task73 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task73(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I2)
   : ta0_(proj), ta1_(I2) { }
};

class Task74 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task74(std::shared_ptr<TATensor<double,4>> I2, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> I3)
   : ta0_(I2), ta1_(Gamma92), ta2_(I3) { }
};

class Task75 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task75(std::shared_ptr<TATensor<double,4>> I3, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I3), ta1_(t2), ta2_(f1) { }
};

class Task76 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task76(std::shared_ptr<TATensor<double,4>> I2, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I6)
   : ta0_(I2), ta1_(f1), ta2_(I6) { }
};

class Task77 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task77(std::shared_ptr<TATensor<double,4>> I6, std::shared_ptr<TATensor<double,6>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I6), ta1_(Gamma2), ta2_(t2) { }
};

class Task78 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task78(std::shared_ptr<TATensor<double,4>> I2, std::shared_ptr<TATensor<double,4>> Gamma3, std::shared_ptr<TATensor<double,4>> I9)
   : ta0_(I2), ta1_(Gamma3), ta2_(I9) { }
};

class Task79 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task79(std::shared_ptr<TATensor<double,4>> I9, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I9), ta1_(t2), ta2_(f1) { }
};

class Task80 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task80(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I11)
   : ta0_(proj), ta1_(I11) { }
};

class Task81 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task81(std::shared_ptr<TATensor<double,4>> I11, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> I12)
   : ta0_(I11), ta1_(Gamma4), ta2_(I12) { }
};

class Task82 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task82(std::shared_ptr<TATensor<double,4>> I12, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I12), ta1_(t2), ta2_(f1) { }
};

class Task83 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task83(std::shared_ptr<TATensor<double,4>> I11, std::shared_ptr<TATensor<double,6>> Gamma5, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I11), ta1_(Gamma5), ta2_(t2) { }
};

class Task84 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task84(std::shared_ptr<TATensor<double,4>> I11, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> I17)
   : ta0_(I11), ta1_(Gamma6), ta2_(I17) { }
};

class Task85 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task85(std::shared_ptr<TATensor<double,4>> I17, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I17), ta1_(t2), e0_(e) { }
};

class Task86 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task86(std::shared_ptr<TATensor<double,4>> I17, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I17), ta1_(t2), ta2_(f1) { }
};

class Task87 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task87(std::shared_ptr<TATensor<double,4>> I17, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I17), ta1_(t2), ta2_(f1) { }
};

class Task88 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task88(std::shared_ptr<TATensor<double,4>> I11, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,2>> I20)
   : ta0_(I11), ta1_(Gamma7), ta2_(I20) { }
};

class Task89 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task89(std::shared_ptr<TATensor<double,2>> I20, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I20), ta1_(t2), ta2_(f1) { }
};

class Task90 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task90(std::shared_ptr<TATensor<double,2>> I20, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I20), ta1_(t2), ta2_(f1) { }
};

class Task91 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task91(std::shared_ptr<TATensor<double,4>> I11, std::shared_ptr<TATensor<double,6>> Gamma9, std::shared_ptr<TATensor<double,4>> I26)
   : ta0_(I11), ta1_(Gamma9), ta2_(I26) { }
};

class Task92 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task92(std::shared_ptr<TATensor<double,4>> I26, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I26), ta1_(t2), ta2_(f1) { }
};

class Task93 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task93(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I31)
   : ta0_(proj), ta1_(I31) { }
};

class Task94 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task94(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> I32)
   : ta0_(I31), ta1_(f1), ta2_(I32) { }
};

class Task95 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task95(std::shared_ptr<TATensor<double,4>> I32, std::shared_ptr<TATensor<double,4>> Gamma3, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I32), ta1_(Gamma3), ta2_(t2) { }
};

class Task96 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task96(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I35)
   : ta0_(I31), ta1_(f1), ta2_(I35) { }
};

class Task97 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task97(std::shared_ptr<TATensor<double,2>> I35, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I35), ta1_(Gamma12), ta2_(t2) { }
};

class Task98 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task98(std::shared_ptr<TATensor<double,4>> I31, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,2>> I38)
   : ta0_(I31), ta1_(f1), ta2_(I38) { }
};

class Task99 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task99(std::shared_ptr<TATensor<double,2>> I38, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I38), ta1_(Gamma12), ta2_(t2) { }
};


}
}
}
#endif
#endif

