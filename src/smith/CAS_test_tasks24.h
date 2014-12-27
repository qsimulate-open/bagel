//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks24.h
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

#ifndef __SRC_SMITH_CAS_test_TASKS24_H
#define __SRC_SMITH_CAS_test_TASKS24_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

class Task1150 : public Task {
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
    Task1150(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1150() {}
};

class Task1151 : public Task {
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
    Task1151(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1151() {}
};

class Task1152 : public Task {
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
    Task1152(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1152() {}
};

class Task1153 : public Task {
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
    Task1153(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1153() {}
};

class Task1154 : public Task {
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
    Task1154(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1154() {}
};

class Task1155 : public Task {
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
    Task1155(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1155() {}
};

class Task1156 : public Task {
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
    Task1156(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1156() {}
};

class Task1157 : public Task {
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
    Task1157(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1157() {}
};

class Task1158 : public Task {
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
    Task1158(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1158() {}
};

class Task1159 : public Task {
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
    Task1159(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1159() {}
};

class Task1160 : public Task {
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
    Task1160(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1160() {}
};

class Task1161 : public Task {
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
    Task1161(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1161() {}
};

class Task1162 : public Task {
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
    Task1162(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1162() {}
};

class Task1163 : public Task {
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
    Task1163(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1163() {}
};

class Task1164 : public Task {
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
    Task1164(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1164() {}
};

class Task1165 : public Task {
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
    Task1165(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1165() {}
};

class Task1166 : public Task {
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
    Task1166(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1166() {}
};

class Task1167 : public Task {
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
    Task1167(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1167() {}
};

class Task1168 : public Task {
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
    Task1168(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1168() {}
};

class Task1169 : public Task {
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
    Task1169(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1169() {}
};

class Task1170 : public Task {
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
    Task1170(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1170() {}
};

class Task1171 : public Task {
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
    Task1171(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1171() {}
};

class Task1172 : public Task {
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
    Task1172(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1172() {}
};

class Task1173 : public Task {
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
    Task1173(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1173() {}
};

class Task1174 : public Task {
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
    Task1174(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1174() {}
};

class Task1175 : public Task {
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
    Task1175(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1175() {}
};

class Task1176 : public Task {
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
    Task1176(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1176() {}
};

class Task1177 : public Task {
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
    Task1177(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1177() {}
};

class Task1178 : public Task {
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
    Task1178(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1178() {}
};

class Task1179 : public Task {
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
    Task1179(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1179() {}
};

class Task1180 : public Task {
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
    Task1180(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1180() {}
};

class Task1181 : public Task {
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
    Task1181(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1181() {}
};

class Task1182 : public Task {
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
    Task1182(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1182() {}
};

class Task1183 : public Task {
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
    Task1183(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1183() {}
};

class Task1184 : public Task {
  protected:
    class Task_local : public SubTask<3,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double e0_;
      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<3,1>(block, in, out), range_(ran), e0_(e) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1184(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e);
    ~Task1184() {}
};

class Task1185 : public Task {
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
    Task1185(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1185() {}
};

class Task1186 : public Task {
  protected:
    class Task_local : public SubTask<3,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
        double e0_;
      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<3,1>(block, in, out), range_(ran), e0_(e) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }
  public:
    Task1186(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e);
    ~Task1186() {}
};

class Task1187 : public Task {
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
    Task1187(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1187() {}
};

class Task1188 : public Task {
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
    Task1188(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1188() {}
};

class Task1189 : public Task {
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
    Task1189(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1189() {}
};

class Task1190 : public Task {
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
    Task1190(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1190() {}
};

class Task1191 : public Task {
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
    Task1191(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1191() {}
};

class Task1192 : public Task {
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
    Task1192(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1192() {}
};

class Task1193 : public Task {
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
    Task1193(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1193() {}
};

class Task1194 : public Task {
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
    Task1194(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1194() {}
};

class Task1195 : public Task {
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
    Task1195(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1195() {}
};

class Task1196 : public Task {
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
    Task1196(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1196() {}
};

class Task1197 : public Task {
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
    Task1197(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1197() {}
};

class Task1198 : public Task {
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
    Task1198(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1198() {}
};

class Task1199 : public Task {
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
    Task1199(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range);
    ~Task1199() {}
};


}
}
}
#endif

