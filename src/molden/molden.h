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

class Molden {
   protected:
      /* This parameter only affects returned values. The functions check the file itself for its formatting. */
      bool is_spherical_;

   private:
      double denormalize(int l, double alpha);

   public:
      Molden(bool is_spherical = true) : is_spherical_(is_spherical) {};
      
      /* Read functions */
      std::vector<std::shared_ptr<Atom> > read_geo(std::string in_file);
      std::shared_ptr<const Coeff> read_mos(std::shared_ptr<const Geometry>, std::string in_file);

      /* Write functions */
      void write_geo(std::shared_ptr<Geometry>, std::string out_file);
      void write_mos(std::shared_ptr<const Reference>, std::string out_file);

      void set_spherical(bool is_spherical) { is_spherical_ = is_spherical; }
};

#endif
