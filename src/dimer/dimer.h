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

#ifndef __dimer_dimer_h
#define __dimer_dimer_h

#include <array>
#include <src/scf/geometry.h>
#include <src/scf/coeff.h>
#include <src/scf/matrix1e.h>
#include <src/wfn/reference.h>

namespace bagel {

typedef std::shared_ptr<const Geometry> RefGeometry;
typedef std::shared_ptr<const Reference> RefReference;
typedef std::shared_ptr<const Coeff> RefCoeff;

/************************************************************************************
*  This class describes a homodimer.                                                *
*************************************************************************************
*  TODO: Extra functionality that would be nice later on:                           *
*     - compute overlap from separate molecules                                     *
************************************************************************************/

class Dimer {
   protected:
      std::pair<RefGeometry,RefGeometry> geompair_;
      std::vector<RefCoeff> coeffs_; // There should only be one or two in this.
      std::shared_ptr<Geometry> supergeometry_;
      std::shared_ptr<Reference> superreference_;
      std::shared_ptr<Coeff> supercoeff_;
      int nbasis_; // Basis size of both together
      
   public:
      Dimer(RefGeometry a, RefGeometry b);
      Dimer(RefGeometry a, std::array<double,3> displacement);
      Dimer(RefReference a, std::array<double,3> displacement);

      RefGeometry get_A() { return geompair_.first; } ;
      RefGeometry get_B() { return geompair_.second; } ;
      std::vector<RefCoeff> coeffs() { return coeffs_; }

      std::shared_ptr<Geometry> supergeom() { return supergeometry_; }
      std::shared_ptr<Reference> superref() { return superreference_; }
      std::shared_ptr<Coeff> supercoeff() { return supercoeff_; }

      int coeffsize() { return coeffs_.size(); }
      int nbasis() {return nbasis_; }

      std::shared_ptr<Coeff> overlap(); 
      void orthonormalize();

      double energy(); // A rather naive, probably temporary, function for computing the energy

   private:
      void construct_geometry();
      void construct_coeff();
};

}

#endif
