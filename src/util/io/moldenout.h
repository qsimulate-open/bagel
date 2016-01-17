//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moldenout.h
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

#ifndef __io_moldenout_h
#define __io_moldenout_h

#include <src/util/io/moldenio.h>

namespace bagel {

class MoldenOut : public MoldenIO {
   protected:
      std::ofstream ofs_;

      void write_geom();
      void write_mos();

   public:
      MoldenOut(std::string filename);

      MoldenOut& operator<< (std::shared_ptr<const Molecule>);
      MoldenOut& operator<< (std::shared_ptr<const Reference>);

};

}

#endif
