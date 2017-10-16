//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gmatrix1e.h
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#ifndef __SRC_MAT1E_GRAD_GMATRIX1E_H
#define __SRC_MAT1E_GRAD_GMATRIX1E_H

#include <src/mat1e/matrix1e.h>

namespace bagel {

template <typename MatType> class GMatrix1eTask;

// specialized matrix for gradient 1e integrals
template <typename MatType = Matrix>
class GMatrix1e {
  friend class GMatrix1eTask<MatType>;
  protected:
    std::vector<std::shared_ptr<MatType>> matrices_;

    virtual void init(std::shared_ptr<const Molecule>);
    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const std::vector<int>&, const std::vector<int>&, std::shared_ptr<const Molecule>) = 0;

  public:
    GMatrix1e() { }
    GMatrix1e(std::shared_ptr<const Molecule>);
    GMatrix1e(const GMatrix1e&);
    virtual ~GMatrix1e() { }

    void ax_plus_y(const double a, const GMatrix1e<MatType>& o) {
      std::transform(o.matrices_.begin(), o.matrices_.end(), matrices_.begin(), matrices_.begin(),
                     [&a](std::shared_ptr<MatType> p, std::shared_ptr<MatType> q) { q->ax_plus_y(a, p); return q; });
    }

    std::shared_ptr<MatType>& data(const int i) { return matrices_[i]; }
    std::shared_ptr<const MatType> data(const int i) const { return matrices_[i]; }

    MatType& operator[](const int i) { return *matrices_[i]; }
    const MatType& operator[](const int i) const { return *matrices_[i]; }

    void fill_upper()       { for (auto& i : matrices_) i->fill_upper(); }
    void fill_upper_conjg() { for (auto& i : matrices_) i->fill_upper_conjg(); }
    template<typename DataType>
    void scale(const DataType a) { for (auto& i : matrices_) i->scale(a); }

    size_t size() const { return matrices_.size(); }

    virtual void print(const std::string name = "", const int len = 0) const;

};

template <typename MatType = Matrix>
class GMatrix1eTask {
  protected:
    GMatrix1e<MatType>* parent_;
    std::array<std::shared_ptr<const Shell>, 2> shell_;
    std::vector<int> atomindex_;
    std::vector<int> offset_;
    std::shared_ptr<const Molecule> mol_;
  public:
    GMatrix1eTask<MatType>(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o, std::shared_ptr<const Molecule> m, GMatrix1e<MatType>* p)
      : parent_(p), shell_(s), atomindex_(a), offset_(o), mol_(m) { }

    void compute() const { parent_->computebatch(shell_, atomindex_, offset_, mol_); }
};

}

extern template class bagel::GMatrix1e<bagel::Matrix>;

#endif
