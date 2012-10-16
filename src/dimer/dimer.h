//
// BAGEL - Parallel electron correlation program.
// Filename: dimer.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __dimer_dimer_h
#define __dimer_dimer_h

#include <array>
#include <src/scf/geometry.h>
#include <src/scf/coeff.h>
#include <src/scf/matrix1e.h>
#include <src/wfn/reference.h>
#include <src/fci/mofile.h>
#include <src/fci/dvec.h>
#include <src/fci/space.h>

namespace bagel {

typedef std::shared_ptr<const Geometry> RefGeometry;
typedef std::shared_ptr<const Reference> RefReference;
typedef std::shared_ptr<const Coeff> RefCoeff;
typedef std::shared_ptr<const Dvec> RefDvec;

/************************************************************************************
*  This class describes a homodimer.                                                *
************************************************************************************/

class Dimer {
   protected:
      std::pair<RefGeometry,RefGeometry> geoms_;
      std::pair<RefReference, RefReference> refs_;

      std::pair<RefCoeff, RefCoeff> coeffs_;
      std::pair<RefDvec, RefDvec> ccvecs_;

      std::shared_ptr<Geometry>   sgeom_;
      std::shared_ptr<Reference>  sref_;
      std::shared_ptr<Coeff>      scoeff_;

      std::shared_ptr<MOFile>     jop_;
      std::pair<std::shared_ptr<MOFile>, std::shared_ptr<MOFile> > jops_; 

      std::shared_ptr<Space>      space_;

      double energy_;      
      std::unique_ptr<double[]> hamiltonian_;

      int dimerbasis_; // Basis size of both together
      int dimerstates_;

      std::pair<int, int> ncore_;
      std::pair<int, int> nact_;
      std::pair<int, int> nvirt_;
      std::pair<int, int> nstates_;
      std::pair<int, int> nbasis_;
      
      const bool symmetric_;
   public:
      // Constructors
      Dimer(RefGeometry a, RefGeometry b);
      Dimer(RefGeometry a, std::array<double,3> displacement);
      Dimer(RefReference a, std::array<double,3> displacement);

      // Return functions
      std::pair<RefGeometry, RefGeometry> geoms() { return geoms_; } ;
      std::pair<RefCoeff, RefCoeff> coeffs() { return coeffs_; } ;
      std::pair<RefDvec, RefDvec> ccvec() { return ccvecs_; } ;

      std::shared_ptr<Geometry> sgeom() { return sgeom_; }
      std::shared_ptr<Reference> sref() { return sref_; }
      std::shared_ptr<Coeff>   scoeff() { return scoeff_; }

      std::pair<const int, const int> nbasis() {return nbasis_; }

      template <int unit> int core(int i) { return (i + unit*ncore_.first); };
      template <int unit> int act(int a)  { return (a + unit*nact_.first); };

      // Utility functions
      std::shared_ptr<Coeff> overlap(); 
      void orthonormalize();

      // Calculations
      void hamiltonian();
      void energy(); // A rather naive, probably temporary, function for computing the energy

   private:
      void construct_geometry();
      void construct_coeff();

      void compute_closeclose();
      void compute_closeactive();
      void compute_activeactive();

      std::pair<std::shared_ptr<Civec>, std::shared_ptr<Civec> > form_sigma_AB(std::pair<std::shared_ptr<Civec>, std::shared_ptr<Civec> > sigma);
      double dot_product(std::pair<std::shared_ptr<Civec>, std::shared_ptr<Civec> > sigma,
                         std::pair<std::shared_ptr<Civec>, std::shared_ptr<Civec> > ccp);
};

}

#endif
