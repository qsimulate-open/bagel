//
// Newint - Parallel electron correlation program.
// Filename: pairfile.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_UTIL_PAIRFILE_H
#define __SRC_UTIL_PAIRFILE_H

#include <tuple>
#include <memory>

template<class T, class U>
class PairFile {
  protected:
    std::shared_ptr<T> file0_;
    std::shared_ptr<U> file1_;

  public:
    // constructors by assignment
    PairFile(std::shared_ptr<T> a, std::shared_ptr<U> b) : file0_(a), file1_(b) {};
    PairFile(std::pair<std::shared_ptr<T>, std::shared_ptr<U> > o) : file0_(o.first), file1_(o.second) {};
    PairFile(std::tuple<std::shared_ptr<T>, std::shared_ptr<U> > o) : file0_(std::get<0>(o)), file1_(std::get<1>(o)) {};
    // copy constructor (that requires T and U to have a copy constructor)
    PairFile(const PairFile<T, U>& o) : file0_(new T(*o.file0_)), file1_(new U(*o.file1_)) {};
    ~PairFile() {};

    std::shared_ptr<T> first() { return file0_; };
    std::shared_ptr<U> second() { return file1_; };
    std::shared_ptr<const T> first() const { return file0_; };
    std::shared_ptr<const U> second() const { return file1_; };

    // operator overloads
    PairFile<T, U> operator+(const PairFile<T, U> o) const { return PairFile<T, U>(*first()+*o->first(), *second()+*o->second()); };
    PairFile<T, U> operator-(const PairFile<T, U> o) const { return PairFile<T, U>(*first()-*o->first(), *second()-*o->second()); };
    PairFile<T, U>& operator+=(const PairFile<T, U> o) { *first()+=*o->first(); *second()+=*o->second(); return *this; };
    PairFile<T, U>& operator-=(const PairFile<T, U> o) { *first()-=*o->first(); *second()-=*o->second(); return *this; };

    // lapack functions
    void daxpy(const double a, const std::shared_ptr<const PairFile<T, U> > o) { first()->daxpy(a, o->first()); second()->daxpy(a, o->second()); }; 
//  double ddot(const std::shared_ptr<PairFile<T, U> > o) const { return first()->ddot(o->first()) + second()->ddot(o->second()); };
    double ddot(const PairFile<T, U>& o) const { return first()->ddot(*o.first()) + second()->ddot(*o.second()); };
    void scale(const double a) { first()->scale(a); second()->scale(a); };

    // assumes that c is already orthogonal with each other.
    double orthog(std::list<std::shared_ptr<const PairFile<T, U> > > c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter) {
        const double scal = - this->ddot(**iter);
        daxpy(scal, **iter);
      }
      const double scal = 1.0/this->norm();
      scale(scal);
      return 1.0/scal; 
    };

};

#endif
