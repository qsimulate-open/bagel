//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimer.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SRC_DIMER_DIMER_H
#define __SRC_DIMER_DIMER_H

#include <src/ci/fci/knowles.h>
#include <src/ci/fci/harrison.h>
#include <src/ci/ras/rasci.h>
#include <src/asd/dimer/dimer_cispace.h>
#include <src/wfn/get_energy.h>
#include <src/util/muffle.h>

namespace bagel {

/// Contains geometries and references for isolated and joined dimers.
/// Used to prepare an ASD calculation and to start the requisite CI calculations

class Dimer : public std::enable_shared_from_this<Dimer> {
  private:
    template <class T>
    using Ref = std::shared_ptr<const T>;

  protected:
    std::shared_ptr<const PTree> input_;

    std::pair<Ref<Geometry>,Ref<Geometry>> geoms_;            ///< Geometry objects of isolated monomers
    std::pair<Ref<Reference>, Ref<Reference>> isolated_refs_; ///< Reference objects of the isolated monomers BEFORE active spaces have been chosen
    std::pair<Ref<Reference>, Ref<Reference>> embedded_refs_; ///< Reference objects for monomers with the other monomer included as an embedding
    std::pair<Ref<Reference>, Ref<Reference>> active_refs_;   ///< Reference objects of the isolated monomers AFTER the active spaces have been chosen

    std::pair<int, int> nvirt_; ///< number of virtual orbitals in sref_ for each unit (needs to be updated when sref_ is

    std::shared_ptr<const Geometry>   sgeom_;                 ///< Supergeometry, i.e., Geometry of dimer
    /// Superreference, i.e., Reference of dimer
    std::shared_ptr<const Reference>  sref_;

    double active_thresh_;                                    ///< overlap threshold for inclusion in the active space
    bool print_orbital_;

  public:
    // Constructors
    Dimer(std::shared_ptr<const PTree> input, Ref<Geometry> a, Ref<Geometry> b); ///< Conjoins the provided Geometry objects
    Dimer(std::shared_ptr<const PTree> input, Ref<Geometry> a); ///< Duplicates provided Geometry according to translation vector specified in input
    Dimer(std::shared_ptr<const PTree> input, Ref<Reference> A, Ref<Reference> B); ///< Conjoins the provided Reference objects
    Dimer(std::shared_ptr<const PTree> input, Ref<Reference> a, const bool linked = false); ///< Duplicates provided Reference according to translation vector specified in input (false) / linked dimer (true)

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

    // Linked dimers
    void set_active(std::shared_ptr<const PTree> idata);
    void reduce_active(std::shared_ptr<const PTree> idata);
    std::shared_ptr<Matrix> form_reference_active_coeff() const;
    std::shared_ptr<Matrix> form_semi_canonical_coeff(std::shared_ptr<const PTree> idata) const;
    std::shared_ptr<Matrix> overlap_selection(std::shared_ptr<const Matrix> control, std::shared_ptr<const Matrix> treatment) const;

    // Calculations
    void scf(std::shared_ptr<const PTree> idata); ///< Driver for preparation of dimer for MultiExcitonHamiltonian or CI calculation

    /// Driver to run monomer CAS calculations used in MEH
    template <class VecType>
    std::shared_ptr<DimerCISpace_base<VecType>> compute_cispace(std::shared_ptr<const PTree> idata);
    /// Driver to run monomer RAS calculations used in MEH
    template <class VecType>
    std::shared_ptr<DimerCISpace_base<VecType>> compute_rcispace(std::shared_ptr<const PTree> idata);

    /// Creates a Reference object for a MEH or ASD calculation
    std::shared_ptr<Reference> build_reference(const int site, const std::vector<bool> meanfield) const;

    // Update
    void update_coeff(std::shared_ptr<const Matrix> coeff);

  private:
    void construct_geometry(const bool linked = false); ///< Forms super geometry (sgeom_) and optionally projects isolated geometries and supergeometry to a specified basis
    void embed_refs();         ///< Forms two references to be used in CI calculations where the inactive monomer is included as "embedding"
    /// Reads information on monomer subspaces from input
    void get_spaces(std::shared_ptr<const PTree> idata, std::vector<std::vector<int>>& spaces_A, std::vector<std::vector<int>>& spaces_B);

    /// Takes monomer references, projections them onto the supergeom basis, and arranges the
    /// to follow (closedA, closedB, activeA, activeB, virtA, virtB) and returns the result
    std::shared_ptr<const Matrix> form_projected_coeffs();

    /// Lowdin orthogonalizes the result of form_projected_coeffs
    std::shared_ptr<const Matrix> construct_coeff(const bool linked = false);

    /// Runs a single set of monomer CAS/RAS calculations as specified
    template <class CIMethod, class VecType>
    std::shared_ptr<const VecType> embedded_ci(std::shared_ptr<const PTree> idata, std::shared_ptr<const Reference> ref,
                                               const int charge, const int nspin, const int nstate, std::string label);
};

template<class CIMethod, class VecType>
std::shared_ptr<const VecType> Dimer::embedded_ci(std::shared_ptr<const PTree> idata, std::shared_ptr<const Reference> ref,
                                              const int charge, const int nspin, const int nstate, std::string label)
{
  auto input = std::make_shared<PTree>(*idata);

  input->erase("charge"); input->put("charge", lexical_cast<std::string>(charge));
  input->erase("nspin"); input->put("nspin", lexical_cast<std::string>(nspin));
  input->erase("ncore"); input->put("ncore", lexical_cast<std::string>(ref->nclosed()));
  input->erase("nstate"); input->put("nstate", lexical_cast<std::string>(nstate));

  // Hiding cout
  std::stringstream outfilename;
  outfilename << "asd_ci_" << label << "_c" << charge << "_s" << nspin;
  Muffle hide(outfilename.str());

  auto ci = std::make_shared<CIMethod>(input, ref->geom(), ref);
  ci->compute();

  return ci->civectors();
}

template <class VecType>
std::shared_ptr<DimerCISpace_base<VecType>> Dimer::compute_cispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  std::pair<int,int> nelea {isolated_refs_.first->nclosed() - active_refs_.first->nclosed(), isolated_refs_.second->nclosed() - active_refs_.second->nclosed()};
  std::pair<int,int> neleb = nelea;

  auto d1 = std::make_shared<Determinants>(active_refs_.first->nact(), nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = std::make_shared<Determinants>(active_refs_.second->nact(), nelea.second, neleb.second, /*compress*/false, /*mute*/true);
  auto out = std::make_shared<DimerCISpace_base<VecType>>(std::make_pair(d1, d2), nelea, neleb);

  std::vector<std::vector<int>> spaces_A, spaces_B;
  get_spaces(idata, spaces_A, spaces_B);

  Timer castime;

  std::shared_ptr<const PTree> fcidata = idata->get_child_optional("fci");
  if (!fcidata) fcidata = std::make_shared<const PTree>();

  auto run_calculations = [this, &fcidata, &castime]
    (std::vector<std::vector<int>> spaces, std::shared_ptr<const Reference> eref, std::shared_ptr<const Reference> aref, std::string label) {
    std::cout << "    Starting embedded CAS-CI calculations on monomer " << label << std::endl;
    std::vector<std::shared_ptr<const VecType>> results;
    for (auto& ispace : spaces) {
      if (ispace.size() != 3) throw std::runtime_error("Spaces should specify \"charge\", \"spin\", and \"nstate\"");

      std::shared_ptr<PTree> input_copy = std::make_shared<PTree>(*fcidata);
      input_copy->erase("norb"); input_copy->put("norb", lexical_cast<std::string>(aref->nact()));

      std::string method = input_copy->get<std::string>("algorithm", "hz");
      std::set<std::string> kh_options = {"kh", "knowles", "handy"};
      std::set<std::string> hz_options = {"hz", "harrison", "zarrabian"};
      std::set<std::string> dist_options = {"dist", "parallel"};

      if (std::find(kh_options.begin(), kh_options.end(), method) != kh_options.end()) {
        using CiType = typename VecType::Ci;
        std::vector<std::shared_ptr<CiType>> tmp;
        std::shared_ptr<const Dvec> vecs = embedded_ci<KnowlesHandy, Dvec>(input_copy, eref, ispace.at(0), ispace.at(1), ispace.at(2), label);
        for (auto& i : vecs->dvec())
          tmp.push_back(std::make_shared<CiType>(*i));
        results.push_back(std::make_shared<VecType>(tmp));
      }
      else if (std::find(hz_options.begin(), hz_options.end(), method) != hz_options.end()) {
        using CiType = typename VecType::Ci;
        std::vector<std::shared_ptr<CiType>> tmp;
        std::shared_ptr<const Dvec> vecs = embedded_ci<HarrisonZarrabian, Dvec>(input_copy, eref, ispace.at(0), ispace.at(1), ispace.at(2), label);
        for (auto& i : vecs->dvec())
          tmp.push_back(std::make_shared<CiType>(*i));
        results.push_back(std::make_shared<VecType>(tmp));
      }
      else
        throw std::runtime_error("Unrecognized FCI type algorithm");

      std::cout << "      - charge: " << ispace.at(0) << ", spin: " << ispace.at(1) << ", nstates: " << ispace.at(2)
                                 << std::fixed << std::setw(10) << std::setprecision(2) << castime.tick() << std::endl;
    }
    return results;
  };

  for (auto& i : run_calculations(spaces_A, embedded_refs_.first, active_refs_.first, "A"))
    out->template insert<0>(i);

  for (auto& i : run_calculations(spaces_B, embedded_refs_.second, active_refs_.second, "B"))
    out->template insert<1>(i);

  return out;
}

template <class VecType>
std::shared_ptr<DimerCISpace_base<VecType>> Dimer::compute_rcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  std::pair<int,int> nelea {isolated_refs_.first->nclosed() - active_refs_.first->nclosed(), isolated_refs_.second->nclosed() - active_refs_.second->nclosed()};
  std::pair<int,int> neleb = nelea;

  // { {nras1, nras2, nras3}, max holes, max particles }
  std::pair<std::tuple<std::array<int, 3>, int, int>, std::tuple<std::array<int, 3>, int, int>> ras_desc;

  // Sample:
  // "restricted" : [ { "orbitals" : [1, 2, 3], "max_holes" : 0, "max_particles" : 2 } ],
  //  puts 1 orbital in RAS1 with no holes allowed, 2 orbital in RAS2, and 3 orbitals in RAS3 with up to 2 particles
  auto restrictions = idata->get_child("restricted");

  auto get_restricted_data = [] (std::shared_ptr<const PTree> i) {
    return std::make_tuple(i->get_array<int, 3>("orbitals"), i->get<int>("max_holes"), i->get<int>("max_particles"));
  };

  if (restrictions->size() == 1) {
    ras_desc = { get_restricted_data(*restrictions->begin()), get_restricted_data(*restrictions->begin()) };
  } else if (restrictions->size() == 2) {
    auto iter = restrictions->begin();
    auto tmp1 = get_restricted_data(*iter++);
    auto tmp2 = get_restricted_data(*iter);
    ras_desc = {tmp1, tmp2};
  } else {
    throw std::logic_error("One or two sets of restrictions must be provided.");
  }

  // This is less than ideal. It'd be better to have some sort of generator object that can be passed around.
  auto d1 = std::make_shared<RASDeterminants>(std::get<0>(ras_desc.first), nelea.first, neleb.first, std::get<1>(ras_desc.first), std::get<2>(ras_desc.first), true);
  auto d2 = std::make_shared<RASDeterminants>(std::get<0>(ras_desc.second), nelea.second, neleb.second, std::get<1>(ras_desc.second), std::get<2>(ras_desc.second), true);

  auto out = std::make_shared<DimerCISpace_base<VecType>>(std::make_pair(d1, d2), nelea, neleb);

  std::vector<std::vector<int>> spaces_A, spaces_B;
  get_spaces(idata, spaces_A, spaces_B);

  Timer rastime;

  std::shared_ptr<const PTree> rasdata = idata->get_child_optional("ras");
  if (!rasdata) rasdata = std::make_shared<const PTree>();

  // Embedded RAS-CI calculations
  auto run_calculations = [this, &rastime, &rasdata] (std::vector<std::vector<int>>& spaces, std::shared_ptr<const Reference> eref,
                                                                                std::tuple<std::array<int, 3>, int, int> ras_desc, std::string label) {
    std::cout << "    Starting embedded RAS-CI calculations on monomer " << label << std::endl;
    std::vector<std::shared_ptr<const VecType>> results;

    for (auto& ispace : spaces) {
      if (ispace.size() != 3)
        throw std::runtime_error("Spaces should specify \"charge\", \"spin\", and \"nstate\"");

      auto input_copy = std::make_shared<PTree>(*rasdata);
      input_copy->erase("max_holes"); input_copy->put("max_holes", lexical_cast<std::string>(std::get<1>(ras_desc)));
      input_copy->erase("max_particles"); input_copy->put("max_particles", lexical_cast<std::string>(std::get<2>(ras_desc)));

      input_copy->erase("active");
      int current = eref->nclosed();
      auto parent = std::make_shared<PTree>();
      for (int i = 0; i < 3; ++i) {
        auto tmp = std::make_shared<PTree>();
        const int nras = std::get<0>(ras_desc)[i];
        for (int i = 0; i < nras; ++i, ++current)
          tmp->push_back(current+1);
        parent->push_back(tmp);
      }
      input_copy->add_child("active", parent);


      std::string method = input_copy->get<std::string>("method", "local");
      std::set<std::string> serial_options = {"local"};
      std::set<std::string> dist_options = {"dist", "parallel"};

      if (std::find(serial_options.begin(), serial_options.end(), method) != serial_options.end())
        results.push_back(std::make_shared<VecType>(embedded_ci<RASCI, RASDvec>(input_copy, eref, ispace.at(0), ispace.at(1), ispace.at(2), label)));
      else
        throw std::logic_error("Unrecognized RAS type algorithm");

      std::cout << "      - charge: " << ispace.at(0) << ", spin: " << ispace.at(1) << ", nstates: " << ispace.at(2)
                                      << std::fixed << std::setw(10) << std::setprecision(2) << rastime.tick() << std::endl;
    }

    return results;
  };

  for (auto& i : run_calculations(spaces_A, embedded_refs_.first, ras_desc.first, "A"))
    out->template insert<0>(i);

  for (auto& i : run_calculations(spaces_B, embedded_refs_.second, ras_desc.second, "B"))
    out->template insert<1>(i);

  return out;
}

}

#endif
