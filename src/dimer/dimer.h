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

#ifndef __SRC_DIMER_DIMER_H
#define __SRC_DIMER_DIMER_H

#include <src/fci/harrison.h>
#include <src/fci/distfci.h>
#include <src/ras/rasci.h>
#include <src/dimer/dimer_cispace.h>
#include <src/wfn/construct_method.h>
#include <src/util/muffle.h>

namespace bagel {

/************************************************************************************
*  This class describes a homodimer.                                                *
************************************************************************************/

class Dimer : public std::enable_shared_from_this<Dimer> {
   template <class T> using Ref = std::shared_ptr<const T>;

   protected:
      std::shared_ptr<const PTree> input_;

      std::pair<Ref<Geometry>,Ref<Geometry>> geoms_;
      std::pair<Ref<Reference>, Ref<Reference>> refs_;
      std::pair<Ref<Reference>, Ref<Reference>> embedded_refs_;
      std::pair<Ref<Coeff>, Ref<Coeff>> coeffs_;

      std::shared_ptr<const Geometry>   sgeom_;
      std::shared_ptr<Reference>  sref_;
      std::shared_ptr<Coeff>      scoeff_;
      std::shared_ptr<Coeff>      proj_coeff_; // Basically the same thing as scoeff_, except purposefully non-orthogonal

      int dimerbasis_; // Basis size of both together
      int nclosed_;

      std::pair<int, int> ncore_;
      std::pair<int, int> nact_;
      std::pair<int, int> nfilledactive_;
      std::pair<int, int> nvirt_;
      std::pair<int, int> nbasis_;
      std::pair<int, int> nele_;

      double active_thresh_;

   public:
      // Constructors
      Dimer(std::shared_ptr<const PTree> input, Ref<Geometry> a, Ref<Geometry> b);
      Dimer(std::shared_ptr<const PTree> input, Ref<Geometry> a);
      Dimer(std::shared_ptr<const PTree> input, Ref<Reference> A, Ref<Reference> B);
      Dimer(std::shared_ptr<const PTree> input, Ref<Reference> a);

      // Return functions
      std::pair<Ref<Geometry>, Ref<Geometry>> geoms() const { return geoms_; };
      std::pair<Ref<Coeff>, Ref<Coeff>> coeffs() const { return coeffs_; };

      std::shared_ptr<const Geometry> sgeom() const { return sgeom_; };
      std::shared_ptr<Reference> sref() const { return sref_; };
      std::shared_ptr<Coeff>   scoeff() const { return scoeff_; };
      std::shared_ptr<Coeff>   proj_coeff() const { return proj_coeff_; };

      void set_sref(std::shared_ptr<const Reference> ref) {
        scoeff_ = std::make_shared<Coeff>(*ref->coeff());
        sref_ = std::make_shared<Reference>(sgeom_, scoeff_, ref->nclosed(), ref->nact(), ref->nvirt());
      }
      void set_coeff(std::shared_ptr<const Matrix> mat) {
        scoeff_ = std::make_shared<Coeff>(*mat);
        sref_->set_coeff(std::make_shared<const Coeff>(*mat));
      };

      std::pair<const int, const int> nbasis() const { return nbasis_; }
      std::pair<const int, const int> ncore() const { return ncore_; }
      std::pair<const int, const int> nact() const { return nact_; }
      std::pair<const int, const int> nfilledactive() const {return nfilledactive_; }
      int dimerbasis() const { return dimerbasis_; }

      // Utility
      void set_active(std::shared_ptr<const PTree> idata);
      void localize(std::shared_ptr<const PTree> idata);

      // Calculations
      void scf(std::shared_ptr<const PTree> idata); // SCF on dimer and then localize
      template <int unit> Ref<Dvec> embedded_casci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates) const;
      template <int unit> Ref<DistDvec> embedded_distcasci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates) const;
      std::shared_ptr<DimerCAS> compute_cispace(std::shared_ptr<const PTree> idata);
      std::shared_ptr<DimerDistCAS> compute_distcispace(std::shared_ptr<const PTree> idata);

      template <int unit> Ref<RASDvec> embedded_rasci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates, std::tuple<std::array<int, 3>, int, int> desc) const;
      std::shared_ptr<DimerRAS> compute_rcispace(std::shared_ptr<const PTree> idata);

   private:
      void construct_geometry();
      void construct_coeff();
      void embed_refs();
};

template<int unit>
std::shared_ptr<const Dvec> Dimer::embedded_casci(const std::shared_ptr<const PTree> idata, const int charge, const int nspin, const int nstate) const {
  const int nclosed = nclosed_;
  const int ncore = (unit == 0) ? nclosed + nfilledactive_.second : nclosed + nfilledactive_.first;
  const int nact = (unit == 0) ? nact_.first : nact_.second;
  const std::shared_ptr<const Reference> embedded_ref = (unit == 0) ? embedded_refs_.first : embedded_refs_.second;

  // Make new input data, set charge, spin to what I want
  auto input = std::make_shared<PTree>(*idata);

  input->erase("charge"); input->put("charge", lexical_cast<std::string>(charge));
  input->erase("nspin"); input->put("nspin", lexical_cast<std::string>(nspin));
  input->erase("ncore"); input->put("ncore", lexical_cast<std::string>(ncore));
  input->erase("norb"); input->put("norb", lexical_cast<std::string>(nact));
  input->erase("nstate"); input->put("nstate", lexical_cast<std::string>(nstate));

  // Hiding normal cout
  std::stringstream outfilename;
  outfilename << "meh_cas_" << (unit == 0 ? "A_c" : "B_c") << charge << "_s" << nspin;
  Muffle hide(outfilename.str());

  std::shared_ptr<FCI> fci;
  fci = std::dynamic_pointer_cast<FCI>(construct_method("fci", input, embedded_ref->geom(), embedded_ref));
  fci->compute();

  return fci->civectors();
}

template<int unit>
std::shared_ptr<const DistDvec> Dimer::embedded_distcasci(const std::shared_ptr<const PTree> idata, const int charge, const int nspin, const int nstate) const {
  const int nclosed = nclosed_;
  const int ncore = (unit == 0) ? nclosed + nfilledactive_.second : nclosed + nfilledactive_.first;
  const int nact = (unit == 0) ? nact_.first : nact_.second;
  const std::shared_ptr<const Reference> embedded_ref = (unit == 0) ? embedded_refs_.first : embedded_refs_.second;

  // Make new input data, set charge, spin to what I want
  auto input = std::make_shared<PTree>(*idata);

  input->erase("charge"); input->put("charge", lexical_cast<std::string>(charge));
  input->erase("nspin"); input->put("nspin", lexical_cast<std::string>(nspin));
  input->erase("ncore"); input->put("ncore", lexical_cast<std::string>(ncore));
  input->erase("norb"); input->put("norb", lexical_cast<std::string>(nact));
  input->erase("nstate"); input->put("nstate", lexical_cast<std::string>(nstate));
  input->erase("algorithm"); input->put("algorithm", "dist");

  // Hiding normal cout
  std::stringstream outfilename;
  outfilename << "meh_distcas_" << (unit == 0 ? "A_c" : "B_c") << charge << "_s" << nspin;
  Muffle hide(outfilename.str());

  std::shared_ptr<DistFCI> fci;
  fci = std::dynamic_pointer_cast<DistFCI>(construct_method("fci", input, embedded_ref->geom(), embedded_ref));
  fci->compute();

  return fci->civectors();
}

template<int unit>
std::shared_ptr<const RASDvec> Dimer::embedded_rasci(const std::shared_ptr<const PTree> idata, const int charge, const int nspin, const int nstate, std::tuple<std::array<int, 3>, int, int> desc) const {
  const int nclosed = nclosed_;
  const int ncore = (unit == 0) ? nclosed + nfilledactive_.second : nclosed + nfilledactive_.first;
  const int nact = (unit == 0) ? nact_.first : nact_.second;
  const std::shared_ptr<const Reference> embedded_ref = (unit == 0) ? embedded_refs_.first : embedded_refs_.second;

  // Make new input data, set charge, spin to what I want
  auto input = std::make_shared<PTree>(*idata);
  auto erase_put = [&input] ( std::string name, int data ) { input->erase(name); input->put(name, lexical_cast<std::string>(data)); };

  erase_put("charge", charge);
  erase_put("nspin", nspin);
  erase_put("ncore", ncore);
  erase_put("nstate", nstate);
  erase_put("max_holes", std::get<1>(desc));
  erase_put("max_particles", std::get<2>(desc));

  input->erase("active");
  int current = ncore;
  auto parent = std::make_shared<PTree>();
  for (int i = 0; i < 3; ++i) {
    auto tmp = std::make_shared<PTree>();
    const int nras = std::get<0>(desc)[i];
    for (int i = 0; i < nras; ++i, ++current)
      tmp->push_back(current+1);
    parent->push_back(tmp);
  }
  input->add_child("active", parent);

  // Hiding normal cout
  std::stringstream outfilename;
  outfilename << "meh_ras_" << (unit == 0 ? "A_c" : "B_c") << charge << "_s" << nspin;
  Muffle hide(outfilename.str());

  std::shared_ptr<RASCI> rasci;
  rasci = std::dynamic_pointer_cast<RASCI>(construct_method("ras", input, embedded_ref->geom(), embedded_ref));
  rasci->compute();

  return rasci->civectors();
}

}

#endif
