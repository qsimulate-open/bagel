
// BAGEL - Parallel electron correlation program.
// Filename: molfile.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __io_moldenio_h
#define __io_moldenio_h

#include <fstream>

#include <src/io/fileio.h>
#include <src/wfn/reference.h>

namespace bagel {

  class MoldenIO : public FileIO {
     protected:
        std::shared_ptr<const Geometry> geom_;
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
