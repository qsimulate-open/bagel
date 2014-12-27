//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks28.h
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

#ifndef __SRC_SMITH_CAS_test_TASKS28_H
#define __SRC_SMITH_CAS_test_TASKS28_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

class Task1350 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1350(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1350() {}
};

class Task1351 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1351(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1351() {}
};

class Task1352 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1352(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1352() {}
};

class Task1353 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1353(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1353() {}
};

class Task1354 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1354(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1354() {}
};

class Task1355 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1355(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1355() {}
};

class Task1356 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1356(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1356() {}
};

class Task1357 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1357(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1357() {}
};

class Task1358 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1358(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1358() {}
};

class Task1359 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1359(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1359() {}
};

class Task1360 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1360(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1360() {}
};

class Task1361 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1361(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1361() {}
};

class Task1362 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1362(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1362() {}
};

class Task1363 : public Task {
  protected:
    class Task_local : public SubTask<1,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1363(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1363() {}
};

class Task1364 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1364(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1364() {}
};

class Task1365 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1365(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1365() {}
};

class Task1366 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1366(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1366() {}
};

class Task1367 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1367(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1367() {}
};

class Task1368 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1368(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1368() {}
};

class Task1369 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1369(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1369() {}
};

class Task1370 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1370(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1370() {}
};

class Task1371 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1371(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1371() {}
};

class Task1372 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1372(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1372() {}
};

class Task1373 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1373(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1373() {}
};

class Task1374 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1374(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1374() {}
};

class Task1375 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1375(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1375() {}
};

class Task1376 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1376(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1376() {}
};

class Task1377 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1377(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1377() {}
};

class Task1378 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1378(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1378() {}
};

class Task1379 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1379(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1379() {}
};

class Task1380 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1380(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1380() {}
};

class Task1381 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1381(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1381() {}
};

class Task1382 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1382(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1382() {}
};

class Task1383 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1383(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1383() {}
};

class Task1384 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1384(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1384() {}
};

class Task1385 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1385(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1385() {}
};

class Task1386 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1386(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1386() {}
};

class Task1387 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1387(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1387() {}
};

class Task1388 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1388(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1388() {}
};

class Task1389 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1389(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1389() {}
};

class Task1390 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1390(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1390() {}
};

class Task1391 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1391(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1391() {}
};

class Task1392 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1392(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1392() {}
};

class Task1393 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1393(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1393() {}
};

class Task1394 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1394(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1394() {}
};

class Task1395 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1395(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1395() {}
};

class Task1396 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1396(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1396() {}
};

class Task1397 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1397(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1397() {}
};

class Task1398 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1398(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1398() {}
};

class Task1399 : public Task {
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
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1399(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1399() {}
};


}
}
}
#endif

