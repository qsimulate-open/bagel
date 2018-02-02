//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_tasks12.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_SMITH_MSCASPT2_TASKS12_H
#define __SRC_SMITH_MSCASPT2_TASKS12_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MSCASPT2{

class Task550 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task550(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task550() {}
};

class Task551 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task551(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task551() {}
};

class Task552 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task552(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task552() {}
};

class Task553 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<4,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<4,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task553(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task553() {}
};

class Task554 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task554(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task554() {}
};

class Task555 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<4,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<4,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task555(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task555() {}
};

class Task556 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task556(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task556() {}
};

class Task557 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<4,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<4,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task557(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task557() {}
};

class Task558 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task558(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task558() {}
};

class Task559 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task559(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task559() {}
};

class Task560 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<6,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<6,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task560(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task560() {}
};

class Task561 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task561(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task561() {}
};

class Task562 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<6,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<6,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task562(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task562() {}
};

class Task563 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task563(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task563() {}
};

class Task564 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask<1,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,1>(std::array<const Index, 1>(), in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task564(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task564() {}
};

class Task565 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task565(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task565() {}
};

class Task566 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask<1,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,1>(std::array<const Index, 1>(), in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task566(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task566() {}
};

class Task567 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task567(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task567() {}
};

class Task568 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<2,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<2,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task568(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task568() {}
};

class Task569 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task569(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task569() {}
};

class Task570 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task570(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task570() {}
};

class Task571 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task571(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task571() {}
};

class Task572 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task572(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task572() {}
};

class Task573 : public Task {  // merged with gamma
  protected:
    std::array<std::shared_ptr<Tensor>,5> out_;
    std::array<std::shared_ptr<const Tensor>,1> in_;
    class Task_local : public SubTask_Merged<4,1,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in,
                   const std::array<std::shared_ptr<Tensor>,5>& out, std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask_Merged<4,1,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_[0]->allocated())
        out_[0]->allocate();
      if (!out_[1]->allocated())
        out_[1]->allocate();
      if (!out_[2]->allocated())
        out_[2]->allocate();
      if (!out_[3]->allocated())
        out_[3]->allocate();
      if (!out_[4]->allocated())
        out_[4]->allocate();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task573(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task573() {}
};

class Task574 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task574(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task574() {}
};

class Task575 : public Task {
  protected:
    std::shared_ptr<Tensor> out_;
    std::array<std::shared_ptr<const Tensor>,2> in_;
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }
        std::shared_ptr<Tensor> out() { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      if (!out_->allocated())
        out_->allocate();
      for (auto& i : in_)
        i->init();
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task575(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task575() {}
};


}
}
}
#endif
#endif

