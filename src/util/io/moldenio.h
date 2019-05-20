//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molfile.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __BAGEL_IO_MOLDENIO_H
#define __BAGEL_IO_MOLDENIO_H

#include <fstream>

#include <src/util/io/fileio.h>
#include <src/wfn/reference.h>

namespace bagel {

class MoldenIO : public FileIO {
  protected:
    std::shared_ptr<const Molecule> mol_;
    std::shared_ptr<const Reference> ref_;

    std::vector<std::vector<int>> m2b_cart_;
    std::vector<std::vector<int>> m2b_sph_;
    std::vector<std::vector<int>> b2m_cart_;
    std::vector<std::vector<int>> b2m_sph_;
    std::vector<std::vector<double>> scaling_;

    void const_scales();
    void const_maps();

    double denormalize(const int l, const double alpha);

  public:
    MoldenIO(std::string filename);
};

namespace molden_impl {
// some implementation utilities
struct complex4 {
  std::complex<double> data[4];

  complex4() { }
  complex4(const complex4& o) = default;
  template<typename T> complex4(const T a) { data[0] = a; data[1] = a; data[2] = a; data[3] = a; } 
  
  complex4& operator=(const complex4& a) { data[0] = a.data[0]; data[1] = a.data[1]; data[2] = a.data[2]; data[3] = a.data[3]; return *this; }
  void operator+=(const complex4& a) { data[0] += a.data[0]; data[1] += a.data[1]; data[2] += a.data[2]; data[3] += a.data[3]; }
  void operator-=(const complex4& a) { data[0] -= a.data[0]; data[1] -= a.data[1]; data[2] -= a.data[2]; data[3] -= a.data[3]; }

  template<typename T> complex4 operator=(const T a) { data[0] = a; data[1] = a; data[2] = a; data[3] = a; return *this; }
  template<typename T> void operator*=(const T& a) { data[0] *= a; data[1] *= a; data[2] *= a; data[3] *= a; }
  template<typename T> void operator/=(const T& a) { data[0] /= a; data[1] /= a; data[2] /= a; data[3] /= a; }
  template<typename T> complex4 operator*(const T& a) const { complex4 out(*this); out *= a; return out; }
  template<typename T> complex4 operator/(const T& a) const { complex4 out(*this); out /= a; return out; }
};
static std::ostream& operator<<(std::ostream& os, const complex4& a) {
  os << std::setw(44) << a.data[0] << " " << std::setw(44) << a.data[1] << " " << std::setw(44) << a.data[2] << " " << std::setw(44) << a.data[3];
  return os;
}

}

}
#endif
