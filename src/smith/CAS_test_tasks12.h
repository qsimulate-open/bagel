//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks12.h
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

#ifndef __SRC_SMITH_CAS_test_TASKS12_H
#define __SRC_SMITH_CAS_test_TASKS12_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

class Task550 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task550(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task550() {}
};

class Task551 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task551(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task551() {}
};

class Task552 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task552(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task552() {}
};

class Task553 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task553(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task553() {}
};

class Task554 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task554(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task554() {}
};

class Task555 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task555(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task555() {}
};

class Task556 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task556(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task556() {}
};

class Task557 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task557(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task557() {}
};

class Task558 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task558(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task558() {}
};

class Task559 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task559(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task559() {}
};

class Task560 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task560(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task560() {}
};

class Task561 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task561(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task561() {}
};

class Task562 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task562(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task562() {}
};

class Task563 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task563(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task563() {}
};

class Task564 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task564(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task564() {}
};

class Task565 : public AccTask {
  protected:
    class Task_local : public SubTask<6,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task565(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task565() {}
};

class Task566 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task566(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task566() {}
};

class Task567 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task567(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task567() {}
};

class Task568 : public AccTask {
  protected:
    class Task_local : public SubTask<2,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task568(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task568() {}
};

class Task569 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task569(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task569() {}
};

class Task570 : public AccTask {
  protected:
    class Task_local : public SubTask<2,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task570(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task570() {}
};

class Task571 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task571(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task571() {}
};

class Task572 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task572(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task572() {}
};

class Task573 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task573(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task573() {}
};

class Task574 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task574(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task574() {}
};

class Task575 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task575(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task575() {}
};

class Task576 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task576(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task576() {}
};

class Task577 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
        double e0_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,1>(block, in, out), range_(ran), e0_(e) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task577(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e);
    ~Task577() {}
};

class Task578 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task578(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task578() {}
};

class Task579 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
        double e0_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,1>(block, in, out), range_(ran), e0_(e) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task579(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e);
    ~Task579() {}
};

class Task580 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task580(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task580() {}
};

class Task581 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task581(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task581() {}
};

class Task582 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task582(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task582() {}
};

class Task583 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task583(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task583() {}
};

class Task584 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task584(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task584() {}
};

class Task585 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task585(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task585() {}
};

class Task586 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task586(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task586() {}
};

class Task587 : public AccTask {
  protected:
    class Task_local : public SubTask<4,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task587(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task587() {}
};

class Task588 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task588(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task588() {}
};

class Task589 : public AccTask {
  protected:
    class Task_local : public SubTask<2,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task589(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task589() {}
};

class Task590 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task590(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task590() {}
};

class Task591 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task591(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task591() {}
};

class Task592 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task592(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task592() {}
};

class Task593 : public AccTask {
  protected:
    class Task_local : public SubTask<6,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task593(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task593() {}
};

class Task594 : public AccTask {
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task594(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task594() {}
};

class Task595 : public AccTask {
  protected:
    class Task_local : public SubTask<6,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task595(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task595() {}
};

class Task596 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task596(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task596() {}
};

class Task597 : public AccTask {
  protected:
    class Task_local : public SubTask<6,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task597(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task597() {}
};

class Task598 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task598(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task598() {}
};

class Task599 : public AccTask {
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double energy_;
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        double energy() const { return energy_; }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      this->target_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->target_ += i->energy();
      }
    }
  public:
    Task599(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task599() {}
};


}
}
}
#endif

