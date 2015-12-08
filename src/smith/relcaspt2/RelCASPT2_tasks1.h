//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_tasks1.h
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

#ifndef __SRC_SMITH_RelCASPT2_TASKS1_H
#define __SRC_SMITH_RelCASPT2_TASKS1_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelCASPT2{

class Task0 : public Task {  // associated with gamma
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
    Task0(std::array<std::shared_ptr<Tensor>,6> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task0() {}
};

class Task1 : public Task {  // associated with gamma
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
    Task1(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task1() {}
};

class Task2 : public Task {  // associated with gamma
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
    Task2(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task2() {}
};

class Task3 : public Task {  // associated with gamma
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
    Task3(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task3() {}
};

class Task4 : public Task {  // associated with gamma
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
    Task4(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task4() {}
};

class Task5 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<8,5> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,8>& block, const std::array<std::shared_ptr<const Tensor>,5>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<8,5>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task5(std::array<std::shared_ptr<Tensor>,6> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task5() {}
};

class Task6 : public Task {  // associated with gamma
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
    Task6(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task6() {}
};

class Task7 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task7(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task7() {}
};

class Task8 : public Task {  // associated with gamma
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
    Task8(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task8() {}
};

class Task9 : public Task {  // associated with gamma
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
    Task9(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task9() {}
};

class Task10 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task10(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task10() {}
};

class Task11 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task11(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task11() {}
};

class Task12 : public Task {  // associated with gamma
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
    Task12(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task12() {}
};

class Task13 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task13(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task13() {}
};

class Task14 : public Task {  // associated with gamma
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
    Task14(std::array<std::shared_ptr<Tensor>,4> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task14() {}
};

class Task15 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task15(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task15() {}
};

class Task16 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task16(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task16() {}
};

class Task17 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task17(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task17() {}
};

class Task18 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,4> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,4>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,4>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task18(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task18() {}
};

class Task19 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task19(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task19() {}
};

class Task20 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task20(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task20() {}
};

class Task21 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<2,1> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor>,1>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task21(std::array<std::shared_ptr<Tensor>,2> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task21() {}
};

class Task22 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task22(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task22() {}
};

class Task23 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task23(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task23() {}
};

class Task24 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task24(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task24() {}
};

class Task25 : public Task {  // associated with gamma
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
    Task25(std::array<std::shared_ptr<Tensor>,5> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task25() {}
};

class Task26 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task26(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task26() {}
};

class Task27 : public Task {  // associated with gamma
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> out_;
    std::shared_ptr<const TATensor<std::complex<double>,4>> rdm2_;
    void compute_() override;
  public:
    Task27(std::shared_ptr<TATensor<std::complex<double>,4>> out, std::shared_ptr<const TATensor<std::complex<double>,4>> rdm2) : out_(out), rdm2_(rdm2) { }
};

class Task28 : public Task {  // associated with gamma
  protected:
    std::shared_ptr<TATensor<std::complex<double>,0>> out_;
    std::shared_ptr<const TATensor<std::complex<double>,2>> rdm1_;
    std::shared_ptr<const TATensor<std::complex<double>,2>> f1_;
    void compute_() override;
  public:
    Task28(std::shared_ptr<TATensor<std::complex<double>,0>> out, std::shared_ptr<const TATensor<std::complex<double>,2>> rdm1,
           std::shared_ptr<const TATensor<std::complex<double>,2>> f1)
      : out_(out), rdm1_(rdm1), f1_(f1) { }
};

class Task29 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task29(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task29() {}
};

class Task30 : public Task {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;
        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }
      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor>,2>& in, std::shared_ptr<Tensor>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2>(block, in, out), range_(ran) { }
        void compute() override;
    };
    std::vector<std::shared_ptr<Task_local>> subtasks_;
    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }
  public:
    Task30(std::array<std::shared_ptr<Tensor>,3> t,  std::array<std::shared_ptr<const IndexRange>,3> range);
    ~Task30() {}
};

class Task31 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task31(std::shared_ptr<TATensor<std::complex<double>,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task32 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task32(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I0)
   : ta0_(proj), ta1_(I0) { }
};

class Task33 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task33(std::shared_ptr<TATensor<std::complex<double>,4>> I0, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I0), ta1_(Gamma0), ta2_(t2) { }
};

class Task34 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task34(std::shared_ptr<TATensor<std::complex<double>,4>> I0, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma94, std::shared_ptr<TATensor<std::complex<double>,4>> I277)
   : ta0_(I0), ta1_(Gamma94), ta2_(I277) { }
};

class Task35 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task35(std::shared_ptr<TATensor<std::complex<double>,4>> I277, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I277), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task36 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task36(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I2)
   : ta0_(proj), ta1_(I2) { }
};

class Task37 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task37(std::shared_ptr<TATensor<std::complex<double>,4>> I2, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma94, std::shared_ptr<TATensor<std::complex<double>,4>> I3)
   : ta0_(I2), ta1_(Gamma94), ta2_(I3) { }
};

class Task38 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task38(std::shared_ptr<TATensor<std::complex<double>,4>> I3, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I3), ta1_(t2), ta2_(f1) { }
};

class Task39 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task39(std::shared_ptr<TATensor<std::complex<double>,4>> I2, std::shared_ptr<TATensor<std::complex<double>,2>> f1, std::shared_ptr<TATensor<std::complex<double>,4>> I6)
   : ta0_(I2), ta1_(f1), ta2_(I6) { }
};

class Task40 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task40(std::shared_ptr<TATensor<std::complex<double>,4>> I6, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I6), ta1_(Gamma2), ta2_(t2) { }
};

class Task41 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task41(std::shared_ptr<TATensor<std::complex<double>,4>> I2, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma3, std::shared_ptr<TATensor<std::complex<double>,4>> I9)
   : ta0_(I2), ta1_(Gamma3), ta2_(I9) { }
};

class Task42 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task42(std::shared_ptr<TATensor<std::complex<double>,4>> I9, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I9), ta1_(t2), ta2_(f1) { }
};

class Task43 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task43(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I11)
   : ta0_(proj), ta1_(I11) { }
};

class Task44 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task44(std::shared_ptr<TATensor<std::complex<double>,4>> I11, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> I12)
   : ta0_(I11), ta1_(Gamma4), ta2_(I12) { }
};

class Task45 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task45(std::shared_ptr<TATensor<std::complex<double>,4>> I12, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I12), ta1_(t2), ta2_(f1) { }
};

class Task46 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task46(std::shared_ptr<TATensor<std::complex<double>,4>> I11, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I11), ta1_(Gamma5), ta2_(t2) { }
};

class Task47 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task47(std::shared_ptr<TATensor<std::complex<double>,4>> I11, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma6, std::shared_ptr<TATensor<std::complex<double>,4>> I17)
   : ta0_(I11), ta1_(Gamma6), ta2_(I17) { }
};

class Task48 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    const double e0_;
    void compute_() override;
  public:
    Task48(std::shared_ptr<TATensor<std::complex<double>,4>> I17, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2, const double e)
   : ta0_(I17), ta1_(t2), ta2_(v2), e0_(e) { }
};

class Task49 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task49(std::shared_ptr<TATensor<std::complex<double>,4>> I17, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> f1)
   : ta0_(I17), ta1_(t2), ta2_(f1) { }
};


}
}
}
#endif
#endif

