//
// Newint - Parallel electron correlation program.
// Filename: quatern.h
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


#ifndef __SRC_UTIL_QUATERN_H
#define __SRC_UTIL_QUATERN_H

#include <array>

// Implements Quaternion for a Dirac equation solver
// T SHOULD NOT BE A POINTER!!

template<typename T>
class Quatern {
  protected:
    std::array<T, 4> data_;

  public:
    Quatern(std::array<T, 4> i) : data_(i) {};
    ~Quatern() {};

    // overloaded operators. Note that T should have these operators

    // addition and subtraction
    Quatern<T> operator+(const Quatern<T>& o) const { return Quatern<T>(std::array<T,4>{{data_[0]+o[0], data_[1]+o[1], data_[2]+o[2], data_[3]+o[3]}}); };
    Quatern<T> operator-(const Quatern<T>& o) const { return Quatern<T>(std::array<T,4>{{data_[0]-o[0], data_[1]-o[1], data_[2]-o[2], data_[3]-o[3]}}); };
    Quatern<T>& operator+=(const Quatern<T>& o) { data_[0] += o[0]; data_[1] += o[1]; data_[2] += o[2]; data_[3] += o[3]; return *this; };
    Quatern<T>& operator-=(const Quatern<T>& o) { data_[0] -= o[0]; data_[1] -= o[1]; data_[2] -= o[2]; data_[3] -= o[3]; return *this; };

    // multiplication
    // (1+i+j+k) * (1'+i'+j'+k) = (11'-ii'-jj'-kk')_1 + (1i'+i1'+jk'-kj')_i + (1j'+j1'+ki'-ik')_j + (1k'+k1'+ij'-ji')_k
    Quatern<T> operator*(const Quatern<T>& o) const { return Quatern<T>(std::array<T,4> {{
      data_[0]*o[0] - data_[1]*o[1] - data_[2]*o[2] - data_[3]*o[3],
      data_[0]*o[1] + data_[1]*o[0] + data_[2]*o[3] - data_[3]*o[2],
      data_[0]*o[2] + data_[2]*o[0] + data_[3]*o[1] - data_[1]*o[3],
      data_[0]*o[3] + data_[3]*o[0] + data_[1]*o[2] - data_[2]*o[1] }}); };
    Quatern<T>& operator*=(const Quatern<T>& o) { Quatern<T> a(*this*o); *this = a; }; 

    T operator[](const size_t i) const { return data_[i]; };

    const T& data(int i) const { return data_[i]; }; 
};

#endif
