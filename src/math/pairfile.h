//
// BAGEL - Parallel electron correlation program.
// Filename: pairfile.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_MATH_PAIRFILE_H
#define __SRC_MATH_PAIRFILE_H

#include <tuple>
#include <cmath>
#include <memory>
#include <list>
#include <src/math/algo.h>

template<class T, class U>
class PairFile {
  protected:
    std::shared_ptr<T> file0_;
    std::shared_ptr<U> file1_;

  public:
    // constructors by assignment
    PairFile(std::shared_ptr<T> a, std::shared_ptr<U> b) : file0_(a), file1_(b) {}
    PairFile(std::pair<std::shared_ptr<T>, std::shared_ptr<U>> o) : file0_(o.first), file1_(o.second) {}
    PairFile(std::tuple<std::shared_ptr<T>, std::shared_ptr<U>> o) : file0_(std::get<0>(o)), file1_(std::get<1>(o)) {}
    // copy constructor (that requires T and U to have a copy constructor)
    PairFile(const PairFile<T, U>& o) : file0_(o.file0_->copy()), file1_(o.file1_->copy()) {}
    ~PairFile() {}

    std::shared_ptr<T>& first() { return file0_; }
    std::shared_ptr<U>& second() { return file1_; }

    std::shared_ptr<const T> first() const { return file0_; }
    std::shared_ptr<const U> second() const { return file1_; }

    // operator overloads
    PairFile<T, U> operator+(const PairFile<T, U>& o) const {
      std::shared_ptr<T> a0 = first()->copy();  *a0 += *o.first();
      std::shared_ptr<U> a1 = second()->copy(); *a1 += *o.second();
      return PairFile(a0, a1);
    }
    PairFile<T, U> operator-(const PairFile<T, U>& o) const {
      std::shared_ptr<T> a0 = first()->copy();  *a0 -= *o.first();
      std::shared_ptr<U> a1 = second()->copy(); *a1 -= *o.second();
      return PairFile(a0, a1);
    }
    PairFile<T, U>& operator+=(const PairFile<T, U>& o) { *first()+=*o.first(); *second()+=*o.second(); return *this; }
    PairFile<T, U>& operator-=(const PairFile<T, U>& o) { *first()-=*o.first(); *second()-=*o.second(); return *this; }

    PairFile<T, U> operator/(const PairFile<T, U>& o) const {
      std::shared_ptr<T> a0 = first()->copy();  *a0 /= *o.first();
      std::shared_ptr<U> a1 = second()->copy(); *a1 /= *o.second();
      return PairFile(a0, a1);
    }
    PairFile<T, U>& operator/=(const PairFile<T, U>& o) { *first()/=*o.first(); *second()/=*o.second(); return *this; }

    // lapack functions
    void ax_plus_y(const double a, const std::shared_ptr<const PairFile<T, U>> o) { first()->ax_plus_y(a, o->first()); second()->ax_plus_y(a, o->second()); }
    double dot_product(const PairFile<T, U>& o) const { return first()->dot_product(*o.first()) + second()->dot_product(*o.second()); }
    double dot_product(const std::shared_ptr<const PairFile<T, U>> o) const { return dot_product(*o); }
    double norm() const { return std::sqrt(dot_product(*this)); }
    double rms() const { return std::sqrt(dot_product(*this)/size()); }
    size_t size() const { return file0_->size() + file1_->size(); }
    void scale(const double a) { first()->scale(a); second()->scale(a); }
    void synchronize() { first()->synchronize(); second()->synchronize(); }

    void zero() { first()->zero(); second()->zero(); }

    std::shared_ptr<PairFile<T, U>> clone() const { return std::make_shared<PairFile<T, U>>(file0_->clone(), file1_->clone()); }

    // assumes that c is already orthogonal with each other.
    double orthog(std::list<std::shared_ptr<const PairFile<T, U>>> c) {
      for (auto& i : c) {
        const double scal = - this->dot_product(*i);
        ax_plus_y(scal, i);
      }
      const double scal = 1.0/this->norm();
      scale(scal);
      return 1.0/scal;
    }

};

#endif
