//
// Newint - Parallel electron correlation program.
// Filename: molden.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __molden_molden_h
#define __molden_molden_h

#include <src/scf/atom.h>
#include <src/scf/geometry.h>
#include <src/wfn/reference.h>
#include <src/util/constants.h>

class Molden {
   protected:
      std::string filename;

   public:
      Molden();
      
      /* Read functions */
      std::vector<std::shared_ptr<Atom> > read_geo(std::string in_file);
      std::shared_ptr<Reference> read_ref(std::string in_file); /* This is still to come */

      /* Write functions, TODO */
      void write_geo(std::shared_ptr<Geometry>, std::string out_file);
      void write_ref(std::shared_ptr<Reference>, std::string out_file);
};

#endif
