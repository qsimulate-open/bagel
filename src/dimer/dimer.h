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

#include <src/fci/distfci.h>
#include <src/ras/rasci.h>
#include <src/ras/distrasci.h>
#include <src/dimer/dimer_cispace.h>
#include <src/wfn/construct_method.h>
#include <src/util/muffle.h>

namespace bagel {

/// Contains geometries and references for isolated and joined dimers.
/// Used to prepare an MEH calculation and to start the requisite CI calculations

class Dimer : public std::enable_shared_from_this<Dimer> {
   template <class T> using Ref = std::shared_ptr<const T>;

   protected:
      std::shared_ptr<const PTree> input_;

      std::pair<Ref<Geometry>,Ref<Geometry>> geoms_;            ///< Geometry objects of isolated monomers
      std::pair<Ref<Reference>, Ref<Reference>> isolated_refs_; ///< Reference objects of the isolated monomers BEFORE active spaces have been chosen
      std::pair<Ref<Reference>, Ref<Reference>> embedded_refs_; ///< Reference objects for monomers with the other monomer included as an embedding
      std::pair<Ref<Reference>, Ref<Reference>> active_refs_;   ///< Reference objects of the isolated monomers AFTER the active spaces have been chosen

      std::shared_ptr<const Geometry>   sgeom_;                 ///< Supergeometry, i.e., Geometry of dimer
      /// Superreference, i.e., Reference of dimer
      /// sref_->coeff() is organized such that the MOs run as (closedA, closedB, activeA, activeB, virtA, virtB)
      std::shared_ptr<const Reference>  sref_;

      double active_thresh_;                                    ///< overlap threshold for inclusion in the active space

   public:
      // Constructors
      Dimer(std::shared_ptr<const PTree> input, Ref<Geometry> a, Ref<Geometry> b); ///< Conjoins the provided Geometry objects
      Dimer(std::shared_ptr<const PTree> input, Ref<Geometry> a); ///< Duplicates provided Geometry according to translation vector specified in input
      Dimer(std::shared_ptr<const PTree> input, Ref<Reference> A, Ref<Reference> B); ///< Conjoins the provided Reference objects
      Dimer(std::shared_ptr<const PTree> input, Ref<Reference> a); ///< Duplicates provided Reference according to translation vector specified in input

      // Return functions
      std::pair<Ref<Geometry>, Ref<Geometry>> geoms() const { return geoms_; };
      std::pair<Ref<Reference>, Ref<Reference>> isolated_refs() const { return isolated_refs_; }
      std::pair<Ref<Reference>, Ref<Reference>> embedded_refs() const { return embedded_refs_; }
      std::pair<Ref<Reference>, Ref<Reference>> active_refs() const { return active_refs_; }

      std::shared_ptr<const Geometry> sgeom() const { return sgeom_; };
      std::shared_ptr<const Reference> sref() const { return sref_; };

      // Utility functions
      /// Sets active space of sref_ using overlaps with isolated_ref_ active spaces
      void set_active(std::shared_ptr<const PTree> idata, const bool localize_first);
      /// Localizes active space and uses the given Fock matrix to diagonalize the subspaces
      void localize(std::shared_ptr<const PTree> idata, std::shared_ptr<const Matrix> fock, const bool localize_first);

      // Calculations
      void scf(std::shared_ptr<const PTree> idata); ///< Driver for preparation of dimer for MultiExcitonHamiltonian or CI calculation

      /// Computes nstates Monomer states for <unit> with specified charge and spin using FCI
      template <int unit> std::shared_ptr<const Dvec> embedded_casci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates) const;
      /// Computes nstates Monomer states for <unit> with specified charge and spin using DistFCI
      template <int unit> std::shared_ptr<const DistDvec> embedded_distcasci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates) const;
      /// Driver to prepare DimerCISpace using FCI
      std::shared_ptr<DimerCAS> compute_cispace(std::shared_ptr<const PTree> idata);
      /// Driver to prepare DimerCISpace using DistFCI
      std::shared_ptr<DimerDistCAS> compute_distcispace(std::shared_ptr<const PTree> idata);

      /// Computes nstates Monomer states for <unit> with specified charge and spin using RASCI
      template <int unit> std::shared_ptr<const RASDvec> embedded_rasci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates, std::tuple<std::array<int, 3>, int, int> desc) const;
      /// Computes nstates Monomer states for <unit> with specified charge and spin using DistRASCI
      template <int unit> std::shared_ptr<const DistRASDvec> embedded_distrasci(std::shared_ptr<const PTree> idata, const int charge, const int spin, const int nstates, std::tuple<std::array<int, 3>, int, int> desc) const;
      /// Driver to prepare DimerCISpace using RASCI
      std::shared_ptr<DimerRAS> compute_rcispace(std::shared_ptr<const PTree> idata);
      /// Driver to prepare DimerCISpace using DistRASCI
      std::shared_ptr<DimerDistRAS> compute_distrcispace(std::shared_ptr<const PTree> idata);

   private:
      void construct_geometry(); ///< Forms super geometry (sgeom_) and optionally projects isolated geometries and supergeometry to a specified basis
      void embed_refs();         ///< Forms two references to be used in CI calculations where the inactive monomer is included as "embedding"
      /// Reads information on monomer subspaces from input
      void get_spaces(std::shared_ptr<const PTree> idata, std::vector<std::vector<int>>& spaces_A, std::vector<std::vector<int>>& spaces_B);

      /// Takes monomer references, projections them onto the supergeom basis, and arranges the
      /// to follow (closedA, closedB, activeA, activeB, virtA, virtB) and returns the result
      std::shared_ptr<const Matrix> form_projected_coeffs();

      /// Lowdin orthogonalizes the result of form_projected_coeffs
      std::shared_ptr<const Matrix> construct_coeff();
};

template<int unit>
std::shared_ptr<const Dvec> Dimer::embedded_casci(const std::shared_ptr<const PTree> idata, const int charge, const int nspin, const int nstate) const {
  const int ncore = (unit == 0) ? active_refs_.first->nclosed() + isolated_refs_.second->nclosed() : active_refs_.second->nclosed() + isolated_refs_.first->nclosed();
  const int nact = (unit == 0) ? active_refs_.first->nact() : active_refs_.second->nact();
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
  const int ncore = (unit == 0) ? active_refs_.first->nclosed() + isolated_refs_.second->nclosed() : active_refs_.second->nclosed() + isolated_refs_.first->nclosed();
  const int nact = (unit == 0) ? active_refs_.first->nact() : active_refs_.second->nact();
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
  const int ncore = (unit == 0) ? active_refs_.first->nclosed() + isolated_refs_.second->nclosed() : active_refs_.second->nclosed() + isolated_refs_.first->nclosed();
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

template<int unit>
std::shared_ptr<const DistRASDvec> Dimer::embedded_distrasci(const std::shared_ptr<const PTree> idata, const int charge, const int nspin, const int nstate, std::tuple<std::array<int, 3>, int, int> desc) const {
  const int ncore = (unit == 0) ? active_refs_.first->nclosed() + isolated_refs_.second->nclosed() : active_refs_.second->nclosed() + isolated_refs_.first->nclosed();
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
  input->erase("algorithm"); input->put("algorithm", "parallel");

  // Hiding normal cout
  std::stringstream outfilename;
  outfilename << "meh_ras_" << (unit == 0 ? "A_c" : "B_c") << charge << "_s" << nspin;
  Muffle hide(outfilename.str());

  std::shared_ptr<DistRASCI> rasci;
  rasci = std::dynamic_pointer_cast<DistRASCI>(construct_method("ras", input, embedded_ref->geom(), embedded_ref));
  rasci->compute();

  return rasci->civectors();
}

}

#endif
