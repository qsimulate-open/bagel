//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradtask.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifdef  GRADTASK_INCLUDE
// *** CAUTION *** this file should be included with GRADTASK_INCLUDE defined (in gradeval_base.h)
// owing to the complicated dependency between GradTask and GradEval_base. Could be cleaned up in the future.

#ifndef __SRC_GRAD_GRADTASK_H
#define __SRC_GRAD_GRADTASK_H

/// Base class for gradient tasks
class GradTask {
  protected:
    std::array<int,4> atomindex_;
    std::array<int,4> offset_;
    GradEval_base* ge_;

    void common_init(const std::vector<int>& a, const std::vector<int>& o) {
      assert(a.size() == o.size());
      int k = 0;
      for (auto i = a.begin(), j = o.begin(); i != a.end(); ++i, ++j, ++k) {
        atomindex_[k] = *i;
        offset_[k] = *j;
      }
    }

  public:
    GradTask(const std::vector<int>& a, const std::vector<int>& o, GradEval_base* p) : ge_(p) { common_init(a,o); }
    virtual void compute() = 0;
};


/// 3-index 2-electron gradient integrals
class GradTask3 : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 4> shell_;
    std::shared_ptr<const DFDist> den_;
  public:
    GradTask3(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
              const std::shared_ptr<const DFDist> d, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), den_(d) { }
    void compute();
};


/// 2-index 2-electron gradient integrals
class GradTask2 : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 4> shell_;
    std::shared_ptr<const Matrix> den2_;
  public:
    GradTask2(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
              const std::shared_ptr<const Matrix> d, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), den2_(d) { }
    void compute();
};


/// 2-index 1-electron gradient integrals
class GradTask1 : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 2> shell_;
    std::shared_ptr<const Matrix> den2_;
    std::shared_ptr<const Matrix> den3_;
    std::shared_ptr<const Matrix> eden_;

    std::shared_ptr<GradFile> compute_nai() const;
    // implemented in gradeval_base.h
    template<typename TBatch>
    std::shared_ptr<GradFile> compute_os(std::shared_ptr<const Matrix> den) const;

  public:
    GradTask1(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o,
              const std::shared_ptr<const Matrix> nmat, const std::shared_ptr<const Matrix> kmat, const std::shared_ptr<const Matrix> omat, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), den2_(nmat), den3_(kmat), eden_(omat) { }
    void compute();
};

/// 2-index 1-electron derivative overlap
class GradTask1s : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 2> shell_;
    std::shared_ptr<const Matrix> den2_;
    std::shared_ptr<const Matrix> den3_;
    std::shared_ptr<const Matrix> eden_;

    // implemented in gradeval_base.h
    template<typename TBatch>
    std::shared_ptr<GradFile> compute_os(std::shared_ptr<const Matrix> den) const;

  public:
    GradTask1s(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o,
               const std::shared_ptr<const Matrix> vmat, const std::shared_ptr<const Matrix> kmat, const std::shared_ptr<const Matrix> omat, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), den2_(omat), den3_(kmat), eden_(vmat) { }
    void compute();
};

/// 2-index 1-electron finite-nucleus NAI gradient integrals
class GradTask1f : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 4> shell_;
    std::shared_ptr<const Matrix> den_;
  public:
    GradTask1f(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
               const std::shared_ptr<const Matrix> d, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), den_(d) { }
    void compute();
};


/// 3-index 2-electron relativistic gradient integrals (small components)
class GradTask3r : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 4> shell_;
    std::array<std::shared_ptr<const DFDist>,6> rden3_;

    std::shared_ptr<GradFile> compute_smalleri() const;
  public:
    GradTask3r(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
               const std::array<std::shared_ptr<const DFDist>,6> rmat, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), rden3_(rmat) { }
    void compute();
};


/// 2-index 1-electron relativistic gradient integrals (small components)
class GradTask1r : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 2> shell_;
    std::array<std::shared_ptr<const Matrix>,6> rden_;

    std::shared_ptr<GradFile> compute_smallnai() const;
  public:
    GradTask1r(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o,
               const std::array<std::shared_ptr<const Matrix>,6> rmat, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), rden_(rmat) { }
    void compute();
};


/// 2-index 1-electron relativistic finite-nucleus NAI gradient integrals
class GradTask1rf : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 4> shell_;
    std::array<std::shared_ptr<const Matrix>,6> rden_;

    std::shared_ptr<GradFile> compute_smalleri() const;
  public:
    GradTask1rf(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
                const std::array<std::shared_ptr<const Matrix>,6> d, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), rden_(d) { }
    void compute();
};


/// 1-electron DKH gradient integrals
class GradTask1d : public GradTask {
  private:
    std::array<std::shared_ptr<const Shell>, 2> shell_;
    std::array<std::shared_ptr<const Matrix>, 4> den_;

    std::shared_ptr<GradFile> compute_nai() const;
    std::shared_ptr<GradFile> compute_smallnai() const;
    // implemented in gradeval_base.h
    template<typename TBatch>
    std::shared_ptr<GradFile> compute_os(std::shared_ptr<const Matrix>) const;

  public:
    GradTask1d(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o,
              const std::array<std::shared_ptr<const Matrix>, 4> den, GradEval_base* p)
      : GradTask(a, o, p), shell_(s), den_(den) { }
    void compute();
};


#endif
#endif
