//
// Newint - Parallel electron correlation program.
// Filename: dimer.h
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

#include <src/scf/geometry.h>
#include <src/scf/coeff.h>
#include <src/scf/matrix1e.h>
#include <src/wfn/reference.h>

#ifndef __dimer_dimer_h
#define __dimer_dimer_h

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<const Reference> RefReference;
typedef std::shared_ptr<const Coeff> RefCoeff;

/************************************************************************************
*  This class describes a homodimer. Since it is a homodimer, there is only one     *
*  important coeff_ matrix. I'll figure out how to generalize the rest of the class *
*  later on, but for now it will just be a way to construct a geometry. The rest    *
*  will come eventually.                                                            *
*************************************************************************************
*  TODO: Extra functionality that would be nice later on:                           *
*     - compute overlap from separate molecules                                     *
*     - output geometry of pair as one object                                       *
*     - output reference of pair as one object                                      *
*     - RDM stuff? No idea if I'll need this or not...                              *
************************************************************************************/

class Dimer {
   protected:
      std::pair<RefGeometry,RefGeometry> geompair_;
      std::vector<RefCoeff> coeffs_; // There should only be one or two in this.
      
   public:
      Dimer(RefGeometry a, RefGeometry b);
      Dimer(RefGeometry a, std::tuple<double,double,double> displacement);
      Dimer(RefReference a, std::tuple<double,double,double> displacement);

      RefGeometry get_A() { return geompair_.first; } ;
      RefGeometry get_B() { return geompair_.second; } ;
      std::vector<RefCoeff> coeffs() { return coeffs_; }

      int coeffsize() { return coeffs_.size(); }

      /* Combine the two geometries into one */
      RefGeometry geometry();
      std::shared_ptr<const Geometry> const_geometry();
      std::shared_ptr<Coeff> coefficients();
      std::shared_ptr<Coeff> coefficients(std::shared_ptr<const Geometry>);
      Matrix1e overlap();
};

#endif
