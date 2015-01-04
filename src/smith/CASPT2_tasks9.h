//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks9.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS9_H
#define __SRC_SMITH_CASPT2_TASKS9_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task400 : public AccTask {
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
    Task400(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task400() {}
};

class Task401 : public AccTask {
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
    Task401(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task401() {}
};

class Task402 : public AccTask {
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
    Task402(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task402() {}
};

class Task403 : public AccTask {
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
    Task403(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task403() {}
};

class Task404 : public AccTask {
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
    Task404(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task404() {}
};

class Task405 : public AccTask {
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
    Task405(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task405() {}
};

class Task406 : public AccTask {
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
    Task406(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task406() {}
};

class Task407 : public AccTask {
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
    Task407(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task407() {}
};

class Task408 : public AccTask {
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
    Task408(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task408() {}
};

class Task409 : public AccTask {
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
    Task409(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task409() {}
};

class Task410 : public AccTask {
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
    Task410(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task410() {}
};

class Task411 : public AccTask {
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
    Task411(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task411() {}
};

class Task412 : public AccTask {
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
    Task412(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task412() {}
};

class Task413 : public AccTask {
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
    Task413(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task413() {}
};

class Task414 : public AccTask {
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
    Task414(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task414() {}
};

class Task415 : public AccTask {
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
    Task415(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task415() {}
};

class Task416 : public AccTask {
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
    Task416(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task416() {}
};

class Task417 : public AccTask {
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
    Task417(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task417() {}
};

class Task418 : public AccTask {
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
    Task418(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task418() {}
};

class Task419 : public AccTask {
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
    Task419(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task419() {}
};

class Task420 : public AccTask {
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
    Task420(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task420() {}
};

class Task421 : public AccTask {
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
    Task421(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task421() {}
};

class Task422 : public AccTask {
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
    Task422(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task422() {}
};

class Task423 : public AccTask {
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
    Task423(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task423() {}
};

class Task424 : public AccTask {
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
    Task424(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task424() {}
};

class Task425 : public AccTask {
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
    Task425(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task425() {}
};

class Task426 : public AccTask {
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
    Task426(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task426() {}
};

class Task427 : public AccTask {
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
    Task427(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task427() {}
};

class Task428 : public AccTask {
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
    Task428(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task428() {}
};

class Task429 : public AccTask {
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
    Task429(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task429() {}
};

class Task430 : public AccTask {
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
    Task430(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task430() {}
};

class Task431 : public AccTask {
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
    Task431(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task431() {}
};

class Task432 : public AccTask {
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
    Task432(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task432() {}
};

class Task433 : public AccTask {
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
    Task433(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task433() {}
};

class Task434 : public AccTask {
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
    Task434(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task434() {}
};

class Task435 : public AccTask {
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
    Task435(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e);
    ~Task435() {}
};

class Task436 : public AccTask {
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
    Task436(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task436() {}
};

class Task437 : public AccTask {
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
    Task437(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e);
    ~Task437() {}
};

class Task438 : public AccTask {
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
    Task438(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task438() {}
};

class Task439 : public AccTask {
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
    Task439(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task439() {}
};

class Task440 : public AccTask {
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
    Task440(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task440() {}
};

class Task441 : public AccTask {
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
    Task441(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task441() {}
};

class Task442 : public AccTask {
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
    Task442(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task442() {}
};

class Task443 : public AccTask {
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
    Task443(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task443() {}
};

class Task444 : public AccTask {
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
    Task444(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task444() {}
};

class Task445 : public AccTask {
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
    Task445(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task445() {}
};

class Task446 : public AccTask {
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
    Task446(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task446() {}
};

class Task447 : public AccTask {
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
    Task447(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task447() {}
};

class Task448 : public AccTask {
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
    Task448(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task448() {}
};

class Task449 : public AccTask {
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
    Task449(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task449() {}
};


}
}
}
#endif

