//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks27.h
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

#ifndef __SRC_SMITH_CAS_test_TASKS27_H
#define __SRC_SMITH_CAS_test_TASKS27_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

class Task1300 : public Task {
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
    Task1300(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1300() {}
};

class Task1301 : public Task {
  protected:
    class Task_local : public SubTask<7,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1301(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1301() {}
};

class Task1302 : public Task {
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
    Task1302(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1302() {}
};

class Task1303 : public Task {
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
    Task1303(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1303() {}
};

class Task1304 : public Task {
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
    Task1304(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1304() {}
};

class Task1305 : public Task {
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
    Task1305(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1305() {}
};

class Task1306 : public Task {
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
    Task1306(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1306() {}
};

class Task1307 : public Task {
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
    Task1307(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1307() {}
};

class Task1308 : public Task {
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
    Task1308(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1308() {}
};

class Task1309 : public Task {
  protected:
    class Task_local : public SubTask<7,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1309(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1309() {}
};

class Task1310 : public Task {
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
    Task1310(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1310() {}
};

class Task1311 : public Task {
  protected:
    class Task_local : public SubTask<7,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double e0_;
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<7,1>(block, in, out), range_(ran), e0_(e) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1311(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e);
    ~Task1311() {}
};

class Task1312 : public Task {
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
    Task1312(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1312() {}
};

class Task1313 : public Task {
  protected:
    class Task_local : public SubTask<7,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1313(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1313() {}
};

class Task1314 : public Task {
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
    Task1314(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1314() {}
};

class Task1315 : public Task {
  protected:
    class Task_local : public SubTask<7,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1315(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1315() {}
};

class Task1316 : public Task {
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
    Task1316(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1316() {}
};

class Task1317 : public Task {
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
    Task1317(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1317() {}
};

class Task1318 : public Task {
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
    Task1318(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1318() {}
};

class Task1319 : public Task {
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
    Task1319(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1319() {}
};

class Task1320 : public Task {
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
    Task1320(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1320() {}
};

class Task1321 : public Task {
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
    Task1321(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1321() {}
};

class Task1322 : public Task {
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
    Task1322(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1322() {}
};

class Task1323 : public Task {
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
    Task1323(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1323() {}
};

class Task1324 : public Task {
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
    Task1324(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1324() {}
};

class Task1325 : public Task {
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
    Task1325(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1325() {}
};

class Task1326 : public Task {
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
    Task1326(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1326() {}
};

class Task1327 : public Task {
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
    Task1327(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1327() {}
};

class Task1328 : public Task {
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
    Task1328(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1328() {}
};

class Task1329 : public Task {
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
    Task1329(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1329() {}
};

class Task1330 : public Task {
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
    Task1330(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1330() {}
};

class Task1331 : public Task {
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
    Task1331(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1331() {}
};

class Task1332 : public Task {
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
    Task1332(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1332() {}
};

class Task1333 : public Task {
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
    Task1333(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1333() {}
};

class Task1334 : public Task {
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
    Task1334(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1334() {}
};

class Task1335 : public Task {
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
    Task1335(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1335() {}
};

class Task1336 : public Task {
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
    Task1336(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1336() {}
};

class Task1337 : public Task {
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
    Task1337(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1337() {}
};

class Task1338 : public Task {
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
    Task1338(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1338() {}
};

class Task1339 : public Task {
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
    Task1339(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1339() {}
};

class Task1340 : public Task {
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
    Task1340(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1340() {}
};

class Task1341 : public Task {
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
    Task1341(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1341() {}
};

class Task1342 : public Task {
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
    Task1342(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1342() {}
};

class Task1343 : public Task {
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
    Task1343(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1343() {}
};

class Task1344 : public Task {
  protected:
    class Task_local : public SubTask<1,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1344(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1344() {}
};

class Task1345 : public Task {
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
    Task1345(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1345() {}
};

class Task1346 : public Task {
  protected:
    class Task_local : public SubTask<1,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1346(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1346() {}
};

class Task1347 : public Task {
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
    Task1347(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1347() {}
};

class Task1348 : public Task {
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
    Task1348(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1348() {}
};

class Task1349 : public Task {
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
    Task1349(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1349() {}
};


}
}
}
#endif

