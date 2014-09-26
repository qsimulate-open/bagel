//
// BAGEL - Parallel electron correlation program.
// Filename: transform.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_TRANSFORM_H
#define __SRC_PERIODIC_TRANSFORM_H

#include <complex>
#include <algorithm>
#include <src/periodic/data.h>
#include <src/periodic/kdata.h>

namespace bagel {

class Transform {
  protected:
    int nbasis_;

    // 3-index objects
    std::shared_ptr<Data> data_;    // (g, i, j)
    std::shared_ptr<KData> kdata_;  // (k, i, j)

    std::vector<std::array<double, 3>> gvector_;
    std::vector<std::array<double, 3>> kvector_;
    int num_gvector_, num_kvector_;

  public:
    Transform(const int, std::shared_ptr<Data>, std::shared_ptr<KData>,
              const std::vector<std::array<double, 3>>, const std::vector<std::array<double, 3>>);
    ~Transform() { }

    // Fourier transform
    void ft();
    // Inverse fourier transform
    void ift();

    const int nbasis() { return nbasis_; }
    const std::shared_ptr<const Data> data() const { return data_; }
    const std::shared_ptr<const KData> kdata() const { return kdata_; }

    const std::vector<std::array<double, 3>> gvector() const { return gvector_; }
    const std::array<double, 3> gvector(const int i) const { return gvector_[i]; }
    const std::vector<std::array<double, 3>> kvector() const { return kvector_; }
    const std::array<double, 3> kvector(const int i) const { return kvector_[i]; }

    const int num_gvector() { return num_gvector_; }
    const int num_kvector() { return num_kvector_; }
};

}

#endif
