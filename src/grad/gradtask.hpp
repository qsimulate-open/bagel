//
// BAGEL - Parallel electron correlation program.
// Filename: gradtask.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#endif
#endif
