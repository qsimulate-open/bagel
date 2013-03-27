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

#include <src/wfn/geometry.h>
#include <src/wfn/reference.h>
#include <src/wfn/ciwfn.h>
#include <src/scf/coeff.h>
#include <src/fci/dvec.h>
#include <src/util/matrix.h>
#include <src/dimer/dimer_cispace.h>

namespace bagel {

typedef std::shared_ptr<const Geometry> RefGeometry;
typedef std::shared_ptr<const Reference> RefReference;
typedef std::shared_ptr<const Coeff> RefCoeff;
typedef std::shared_ptr<const Dvec> RefDvec;
typedef std::shared_ptr<const CIWfn> RefCIWfn;
typedef std::multimap<std::string,std::string> MultimapInput;

/************************************************************************************
*  This class describes a homodimer.                                                *
************************************************************************************/

class Dimer : public std::enable_shared_from_this<Dimer> {
   protected:
      std::pair<RefGeometry,RefGeometry> geoms_;
      std::pair<RefReference, RefReference> refs_;
      std::pair<RefReference, RefReference> embedded_refs_;

      std::pair<RefCoeff, RefCoeff> coeffs_;
      std::pair<RefDvec, RefDvec> ccvecs_;
      std::pair<RefCIWfn, RefCIWfn> ci_;

      std::shared_ptr<Geometry>   sgeom_;
      std::shared_ptr<Reference>  sref_;
      std::shared_ptr<Coeff>      scoeff_;
      std::shared_ptr<Coeff>      proj_coeff_; // Basically the same thing as scoeff_, except purposefully non-orthogonal

      int dimerbasis_; // Basis size of both together
      int dimerstates_;

      std::pair<int, int> ncore_;
      std::pair<int, int> nact_;
      std::pair<int, int> nfilledactive_;
      std::pair<int, int> nvirt_;
      std::pair<int, int> nstates_;
      std::pair<int, int> nbasis_;
      std::pair<int, int> nele_;

   public:
      // Constructors
      Dimer(RefGeometry a, RefGeometry b);
      Dimer(RefGeometry a, std::array<double,3> displacement);
      Dimer(RefReference A, RefReference B);
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

      void set_sref(std::shared_ptr<const Reference> ref) { 
        scoeff_ = std::shared_ptr<Coeff>(new Coeff(*ref->coeff()));
        sref_ = std::shared_ptr<Reference>(new Reference(sgeom_, scoeff_, ref->nclosed(), ref->nact(), ref->nvirt()));
      }
      void set_coeff(std::shared_ptr<const Coeff> coeff) {
        scoeff_ = std::shared_ptr<Coeff>(new Coeff(*coeff));
        sref_->set_coeff(coeff);
      };
      void set_coeff(std::shared_ptr<const Matrix> mat) {
        scoeff_ = std::shared_ptr<Coeff>(new Coeff(*mat));
        sref_->set_coeff(std::shared_ptr<const Coeff>(new const Coeff(*mat)));
      };

      std::pair<const int, const int> nbasis() const { return nbasis_; }
      std::pair<const int, const int> ncore() const { return ncore_; }
      std::pair<const int, const int> nact() const { return nact_; }
      std::pair<const int, const int> nfilledactive() const {return nfilledactive_; }
      int dimerbasis() const { return dimerbasis_; }

      // Utility functions
      std::shared_ptr<Coeff> overlap() const; 
      std::shared_ptr<Matrix> form_density_rhf(std::shared_ptr<const Coeff> coeff) const { return std::shared_ptr<Matrix>(); };

      void set_active(MultimapInput idata);
      void localize(MultimapInput idata);


      // Calculations
      void scf(MultimapInput idata); // SCF on dimer and then localize
      std::pair<RefDvec,RefDvec> embedded_casci(MultimapInput idata, const int charge, const int spin, const int nstates) const;
      std::shared_ptr<DimerCISpace> compute_cispace(MultimapInput idata);

   private:
      void construct_geometry();
      void construct_coeff();
      void embed_refs();
      void orthonormalize();
};

}

#endif
