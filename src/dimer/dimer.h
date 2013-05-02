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
#include <src/fci/harrison.h>
#include <src/util/matrix.h>
#include <src/dimer/dimer_cispace.h>

namespace bagel {

/************************************************************************************
*  This class describes a homodimer.                                                *
************************************************************************************/

class Dimer : public std::enable_shared_from_this<Dimer> {
   template <class T> using Ref = std::shared_ptr<const T>;
   using MultimapInput = std::multimap<std::string,std::string>;

   protected:
      std::pair<Ref<Geometry>,Ref<Geometry>> geoms_;
      std::pair<Ref<Reference>, Ref<Reference>> refs_;
      std::pair<Ref<Reference>, Ref<Reference>> embedded_refs_;

      std::pair<Ref<Coeff>, Ref<Coeff>> coeffs_;
      std::pair<Ref<Dvec>, Ref<Dvec>> ccvecs_;
      std::pair<Ref<CIWfn>, Ref<CIWfn>> ci_;

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
      Dimer(Ref<Geometry> a, Ref<Geometry> b);
      Dimer(Ref<Geometry> a, std::array<double,3> displacement);
      Dimer(Ref<Reference> A, Ref<Reference> B);
      Dimer(Ref<Reference> a, std::array<double,3> displacement);
      Dimer(Ref<CIWfn> a, std::array<double,3> displacement);

      // Return functions
      std::pair<Ref<Geometry>, Ref<Geometry>> geoms() const { return geoms_; };
      std::pair<Ref<Coeff>, Ref<Coeff>> coeffs() const { return coeffs_; };
      std::pair<Ref<Dvec>, Ref<Dvec>> ccvec() const { return ccvecs_; };

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
      template <int unit> Ref<Dvec> embedded_casci(MultimapInput idata, const int charge, const int spin, const int nstates) const;
      std::shared_ptr<DimerCISpace> compute_cispace(MultimapInput idata);

   private:
      void construct_geometry();
      void construct_coeff();
      void embed_refs();
      void orthonormalize();
};

template<int unit>
std::shared_ptr<const Dvec> Dimer::embedded_casci(std::multimap<std::string,std::string> idata, const int charge, const int nspin, const int nstates) const { 
  const int nclosed = ncore_.first + ncore_.second;
  const int ncore = (unit == 0) ? nclosed + nfilledactive_.second : nclosed + nfilledactive_.first;
  const int nact = (unit == 0) ? nact_.first : nact_.second;
  const std::shared_ptr<const Reference> embedded_ref = (unit == 0) ? embedded_refs_.first : embedded_refs_.second;

  // Make new input data, set charge and spin to what I want
  std::multimap<std::string,std::string> input = idata;
  input.erase("charge"); input.erase("nspin");
  input.insert(std::make_pair(std::string("charge"), lexical_cast<std::string>(charge)));
  input.insert(std::make_pair(std::string("nspin"), lexical_cast<std::string>(nspin)));

  std::shared_ptr<FCI> fci(new HarrisonZarrabian(input, embedded_ref, ncore, nact, nstates));
  fci->compute();

  return fci->civectors();
}

}

#endif
