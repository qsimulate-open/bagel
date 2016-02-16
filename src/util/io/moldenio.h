//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molfile.h
// Copyright (C) 2012 Shane Parker
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
}
#endif
