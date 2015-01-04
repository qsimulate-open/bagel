//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks13.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS13_H
#define __SRC_SMITH_CASPT2_TASKS13_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task600 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task600(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task600() {}
};

class Task601 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task601(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task601() {}
};

class Task602 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task602(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task602() {}
};

class Task603 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task603(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task603() {}
};

class Task604 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task604(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task604() {}
};

class Task605 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task605(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task605() {}
};

class Task606 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task606(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task606() {}
};

class Task607 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task607(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task607() {}
};

class Task608 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task608(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task608() {}
};

class Task609 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task609(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task609() {}
};

class Task610 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task610(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task610() {}
};

class Task611 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task611(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task611() {}
};

class Task612 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task612(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task612() {}
};

class Task613 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task613(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task613() {}
};

class Task614 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task614(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task614() {}
};

class Task615 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task615(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task615() {}
};

class Task616 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task616(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task616() {}
};

class Task617 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task617(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task617() {}
};

class Task618 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task618(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task618() {}
};

class Task619 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task619(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task619() {}
};

class Task620 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task620(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task620() {}
};

class Task621 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task621(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task621() {}
};

class Task622 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task622(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task622() {}
};

class Task623 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task623(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task623() {}
};

class Task624 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
        double e0_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,1>(block, in, out), range_(ran), e0_(e) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task624(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e);
    ~Task624() {}
};

class Task625 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task625(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task625() {}
};

class Task626 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task626(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task626() {}
};

class Task627 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task627(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task627() {}
};

class Task628 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task628(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task628() {}
};

class Task629 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task629(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task629() {}
};

class Task630 : public AccTask {
  protected:
    class Task_local : public SubTask<2,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task630(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task630() {}
};

class Task631 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task631(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task631() {}
};

class Task632 : public AccTask {
  protected:
    class Task_local : public SubTask<2,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task632(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task632() {}
};

class Task633 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task633(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task633() {}
};

class Task634 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task634(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task634() {}
};

class Task635 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task635(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task635() {}
};

class Task636 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task636(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task636() {}
};

class Task637 : public AccTask {
  protected:
    class Task_local : public SubTask<2,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task637(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task637() {}
};

class Task638 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task638(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task638() {}
};

class Task639 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task639(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task639() {}
};

class Task640 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task640(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task640() {}
};

class Task641 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task641(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task641() {}
};

class Task642 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task642(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task642() {}
};

class Task643 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
        double e0_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,1>(block, in, out), range_(ran), e0_(e) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task643(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e);
    ~Task643() {}
};

class Task644 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task644(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task644() {}
};

class Task645 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task645(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task645() {}
};

class Task646 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task646(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task646() {}
};

class Task647 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task647(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task647() {}
};

class Task648 : public AccTask {
  protected:
    class Task_local : public SubTask<1,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,2>(std::array<const Index, 1>(), in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task648(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task648() {}
};

class Task649 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double target_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double target() const { return target_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->target();
      }
    }
  public:
    Task649(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task649() {}
};


}
}
}
#endif

