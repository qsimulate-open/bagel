
// BAGEL - Parallel electron correlation program.
// Filename: moldenout.h
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

#ifndef __io_moldenout_h
#define __io_moldenout_h

#include <src/io/moldenio.h>

namespace bagel {

class MoldenOut : public MoldenIO {
   protected:
      std::ofstream ofs_;

      void write_geom();
      void write_mos();
      void write_rel_mos();

   public:
      MoldenOut(std::string filename);

      MoldenOut& operator<< (std::shared_ptr<const Molecule>);
      MoldenOut& operator<< (std::shared_ptr<const Reference>);

};

}

#endif
