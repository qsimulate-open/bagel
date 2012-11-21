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
#include <src/wfn/ciwfn.h>
#include <src/dimer/dimer_jop.h>
#include <src/fci/dvec.h>
#include <src/fci/space.h>
#include <src/util/matrix.h>

namespace bagel {

typedef std::shared_ptr<const Geometry> RefGeometry;
typedef std::shared_ptr<const Reference> RefReference;
typedef std::shared_ptr<const Coeff> RefCoeff;
typedef std::shared_ptr<const Dvec> RefDvec;
typedef std::shared_ptr<const CIWfn> RefCIWfn;

/************************************************************************************
*  This class describes a homodimer.                                                *
************************************************************************************/

class Dimer {
   protected:
      std::pair<RefGeometry,RefGeometry> geoms_;
      std::pair<RefReference, RefReference> refs_;

      std::pair<RefCoeff, RefCoeff> coeffs_;
      std::pair<RefDvec, RefDvec> ccvecs_;
      std::pair<RefCIWfn, RefCIWfn> ci_;

      std::shared_ptr<Geometry>   sgeom_;
      std::shared_ptr<Reference>  sref_;
      std::shared_ptr<Coeff>      scoeff_;
      std::shared_ptr<Coeff>      proj_coeff_; // Basically the same thing as scoeff_, except purposefully non-orthogonal

      std::shared_ptr<DimerJop> jop_;

      std::shared_ptr<Space>      space_;

      double energy_;      
      std::shared_ptr<Matrix> hamiltonian_;

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
      Dimer(RefCIWfn a, std::array<double,3> displacement);

      // Return functions
      std::pair<RefGeometry, RefGeometry> geoms() const { return geoms_; };
      std::pair<RefCoeff, RefCoeff> coeffs() const { return coeffs_; };
      std::pair<RefDvec, RefDvec> ccvec() const { return ccvecs_; };

      std::shared_ptr<Geometry> sgeom() const { return sgeom_; };
      std::shared_ptr<Reference> sref() const { return sref_; };
      std::shared_ptr<Coeff>   scoeff() const { return scoeff_; };
      std::shared_ptr<Coeff>   proj_coeff() const { return proj_coeff_; };

      std::pair<const int, const int> nbasis() const {return nbasis_; };
      std::pair<const int, const int> ncore() const { return ncore_; };
      int dimerbasis() const { return dimerbasis_; };

      template <int unit> int core(int i) const { return (i + unit*ncore_.first); };
      template <int unit> int act(int a) const { return (a + unit*nact_.first); };

      // Utility functions
      std::shared_ptr<Coeff> overlap() const; 
      std::shared_ptr<Matrix> form_density_rhf(std::shared_ptr<const Coeff> coeff) const;

      // Calculations
      void hamiltonian();
      void energy(); // A rather naive, probably temporary, function for computing the energy

   private:
      void construct_geometry();
      void construct_coeff();
      void orthonormalize();

      std::shared_ptr<Matrix> compute_closeclose();
      std::shared_ptr<Matrix> compute_closeactive();
      std::shared_ptr<Matrix> compute_intra_activeactive();
      std::shared_ptr<Matrix> compute_inter_activeactive();

      std::shared_ptr<Dvec> form_sigma_1e(std::shared_ptr<const Dvec> ccvec, double* hdata, const int ij) const;
      std::shared_ptr<Dvec> form_sigma_2e(std::shared_ptr<const Dvec> ccvec, double* mo2e_ptr, const int nact) const;

      void sigma_2aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const;
      void sigma_2bb(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const;
      void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d, const int nact) const;
      void sigma_2ab_2(std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, double* mo2e_ptr) const;
      void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e, const int nact) const;
      
      std::shared_ptr<Matrix> form_EFmatrices_alpha(std::shared_ptr<const Dvec> ccvec, const int ij, const int nstates) const;
      std::shared_ptr<Matrix> form_EFmatrices_beta(std::shared_ptr<const Dvec> ccvec, const int ij, const int nstates) const;

      std::shared_ptr<Matrix> form_JKmatrix(const int ijA, const int ijB) const;
      std::shared_ptr<Matrix> form_Jmatrix(const int ijA, const int ijB) const;
};

}

#endif
