//
// BAGEL - Parallel electron correlation program.
// Filename: quatern.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __SRC_UTIL_QUATERN_H
#define __SRC_UTIL_QUATERN_H

#include <cassert>
#include <cmath>
#include <array>
#include <iomanip>
#include <iostream>
#include <initializer_list>

// Implements Quaternion for a Dirac equation solver
// T SHOULD NOT BE A POINTER!!

namespace bagel {

template<typename T>
class Quatern {
  protected:
    std::array<T, 4>__attribute__((aligned(32))) data_;

  public:
    Quatern(std::array<T, 4> i) : data_(i) {};
    Quatern(std::initializer_list<T> list) {
      auto j = data_.begin();
      for (auto i = list.begin(); i != list.end(); ++i, ++j) *j = *i;
    };
    Quatern(std::array<T, 3> i) {
      data_[0] = T();
      data_[1] = i[0];
      data_[2] = i[1];
      data_[3] = i[2];
    };
    ~Quatern() {};

    // overloaded operators. Note that T should have these operators

    // addition and subtraction
    Quatern<T> operator+(const Quatern<T>& o) const { return Quatern<T>{{data_[0]+o[0], data_[1]+o[1], data_[2]+o[2], data_[3]+o[3]}}; };
    Quatern<T> operator-(const Quatern<T>& o) const { return Quatern<T>{{data_[0]-o[0], data_[1]-o[1], data_[2]-o[2], data_[3]-o[3]}}; };
    Quatern<T>& operator+=(const Quatern<T>& o) { data_[0] += o[0]; data_[1] += o[1]; data_[2] += o[2]; data_[3] += o[3]; return *this; };
    Quatern<T>& operator-=(const Quatern<T>& o) { data_[0] -= o[0]; data_[1] -= o[1]; data_[2] -= o[2]; data_[3] -= o[3]; return *this; };

    // multiplication
    // (1+i+j+k) * (1'+i'+j'+k) = (11'-ii'-jj'-kk')_1 + (1i'+i1'+jk'-kj')_i + (1j'+j1'+ki'-ik')_j + (1k'+k1'+ij'-ji')_k
    Quatern<T> operator*(const Quatern<T>& o) const { return Quatern<T>{{
      data_[0]*o[0] - data_[1]*o[1] - data_[2]*o[2] - data_[3]*o[3],
      data_[0]*o[1] + data_[1]*o[0] + data_[2]*o[3] - data_[3]*o[2],
      data_[0]*o[2] + data_[2]*o[0] + data_[3]*o[1] - data_[1]*o[3],
      data_[0]*o[3] + data_[3]*o[0] + data_[1]*o[2] - data_[2]*o[1] }}; };
    Quatern<T>& operator*=(const Quatern<T>& o) { Quatern<T> a(*this*o); *this = a; return *this; };

    Quatern<T>& operator*=(const double& a) { data_[0] *= a; data_[1] *= a; data_[2] *= a; data_[3] *= a; return *this; };
    Quatern<T> operator*(const double& a) const { Quatern<T> out(*this); out *= a; return out; };
    Quatern<T>& operator/=(const double& a) { data_[0] /= a; data_[1] /= a; data_[2] /= a; data_[3] /= a; return *this; };
    Quatern<T> operator/(const double& a) const { Quatern<T> out(*this); out /= a; return out; };

    const T& data(int i) const { return data_[i]; };
    T& operator[](const int i) { return data_[i]; };
    const T& operator[](const int i) const { return data_[i]; };

    Quatern<T> dagger() const {
      return Quatern<T>{{ data_[0], data_[1]*(-1.0), data_[2]*(-1.0), data_[3]*(-1.0) }};
    }

    std::array<T, 3> ijk() { return std::array<T, 3>{{data_[1], data_[2], data_[3]}}; };

    void print() const { };
    double norm() const { assert(false); };
    double ddot(const Quatern<T>&) const { assert(false); };
    void normalize() { assert(false); };
};

template<> void Quatern<double>::print() const;
template<> double Quatern<double>::norm() const;
template<> double Quatern<double>::ddot(const Quatern<double>& o) const;
template<> void Quatern<double>::normalize();

}

#endif
