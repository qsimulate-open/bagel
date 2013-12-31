//
// BAGEL - Parallel electron correlation program.
// Filename: meh.h
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

#ifndef __MEH_MEH_H
#define __MEH_MEH_H

#include <src/dimer/dimer.h>
#include <src/dimer/dimer_jop.h>
#include <src/dimer/dimer_prop.h>
#include <src/math/davidson.h>
#include <src/meh/gamma_forest.h>
#include <src/meh/meh_spin.h>

namespace bagel {

/************************************************************************************
*  This class computes the Hamiltonian matrix for multiexciton states described by  *
* CAS calculations                                                                  *
************************************************************************************/

enum class Coupling {
  none = 0,
  diagonal = 1,
  aET = 2,
  inv_aET = -2,
  bET = 3,
  inv_bET = -3,
  abFlip = 4,
  baFlip = -4,
  abET = 5,
  inv_abET = -5,
  aaET = 6,
  inv_aaET = -6,
  bbET = 7,
  inv_bbET = -7
};

template <class VecType>
class DimerSubspace_base { // until I come up with a better name
  protected:
    const int offset_;
    const int nstatesA_;
    const int nstatesB_;
    const std::string stringA_;
    const std::string stringB_;

    std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> ci_;

  public:
    DimerSubspace_base(int& _offset, const SpaceKey Akey, const SpaceKey Bkey, std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> _ci) :
      offset_(_offset), nstatesA_(_ci.first->ij()), nstatesB_(_ci.second->ij()), stringA_(Akey.to_string()), stringB_(Bkey.to_string()),
       ci_(_ci) { _offset += dimerstates(); }

    const int offset() const { return offset_; }
    const int dimerstates() const { return nstatesA_ * nstatesB_; }
    const int dimerindex(const int iA, const int iB) const { return (iA + iB*nstatesA_); }
    const std::string string(const int i, const int j) const {
      std::string out = stringA_ + lexical_cast<std::string>(i) + std::string(" ") + stringB_ + lexical_cast<std::string>(j);
      return out;
    }

    template <int unit> const int nstates() const { return ( unit == 0 ? nstatesA_ : nstatesB_ ); }
    template <int unit> std::shared_ptr<const VecType> ci() const { return ( unit == 0 ? ci_.first : ci_.second ); }

};

template <class VecType>
class MultiExcitonHamiltonian {
   protected: using DSubSpace = DimerSubspace_base<VecType>;
   protected: using DCISpace = DimerCISpace_base<VecType>;
   protected: using CiType = typename VecType::Ci;

   protected:
      std::shared_ptr<const Dimer> dimer_;
      std::shared_ptr<const Reference> ref_;

      std::shared_ptr<DimerJop> jop_;

      std::shared_ptr<DCISpace> cispace_;
      std::shared_ptr<GammaForest<VecType, 2>> gammaforest_;
      std::vector<DSubSpace> subspaces_;

      std::unique_ptr<double[]> denom_;

      std::shared_ptr<Matrix> hamiltonian_;
      std::shared_ptr<Matrix> adiabats_; // Eigenvectors of adiabatic states
      std::vector<std::pair<std::string, std::shared_ptr<Matrix>>> properties_;

      int max_spin_;
      std::shared_ptr<MEHSpin> spin_;

      std::vector<double> energies_; // Adiabatic energies

      // Total system quantities
      const int dimerbasis_;
      const int dimerclosed_;
      const int dimeractive_;
      int dimerstates_;

      // Localized quantities
      std::pair<const int, const int> nact_;
      std::pair<const int, const int> nbasis_;

      // Options
      int nstates_;
      int nspin_;
      int max_iter_;
      int davidsonceiling_;
      bool store_matrix_;
      bool dipoles_;

      double thresh_;
      double print_thresh_;

   public:
      MultiExcitonHamiltonian(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DCISpace> cispace);

      void compute();

      std::vector<double> energy() const { return energies_; }
      double energy(const int i) const { return energies_.at(i); }

      void print_hamiltonian(const std::string title = "MultiExciton Hamiltonian", const int nstates = 10) const;
      void print_adiabats(const double thresh = 0.05, const std::string title = "Adiabats", const int nstates = 10) const;
      void print_property(const std::string label, std::shared_ptr<const Matrix>, const int size = 10) const ;
      void print(const int nstates = 10, const double thresh = 0.05) const;

      const Coupling coupling_type(const DSubSpace& AB, const DSubSpace& ApBp) const;

   private:
      void reorder_matrix(const double* source, double* target, const int nA, const int nAp, const int nB, const int nBp) const;

      std::shared_ptr<const Matrix> apply_hamiltonian(const Matrix& o);

      std::shared_ptr<Matrix> compute_1e_prop(std::shared_ptr<const Matrix> hAA, std::shared_ptr<const Matrix> hBB, std::shared_ptr<const Matrix> hAB, const double core) const;
      std::shared_ptr<Matrix> compute_offdiagonal_1e(const DSubSpace& AB, const DSubSpace& ApBp, std::shared_ptr<const Matrix> hAB) const;
      std::shared_ptr<Matrix> compute_diagonal_1e(const DSubSpace& subspace, const double* hAA, const double* hBB, const double diag) const;

      // Diagonal block stuff
      std::shared_ptr<Matrix> compute_diagonal_block(DSubSpace& subspace);
      std::shared_ptr<Matrix> compute_intra(const DSubSpace& subspace, std::shared_ptr<const DimerJop> jop, const double diag);

      virtual std::shared_ptr<VecType> form_sigma(std::shared_ptr<const VecType> ccvec, std::shared_ptr<const MOFile> jop) const = 0;
      virtual std::shared_ptr<VecType> form_sigma_1e(std::shared_ptr<const VecType> ccvec, const double* modata) const = 0;

      // Gamma Tree building
      void gamma_couple_blocks(DSubSpace& AB, DSubSpace& ApBp);
      void spin_couple_blocks(DSubSpace& AB, DSubSpace& ApBp, std::map<std::pair<int, int>, double>& spinmap); // Off-diagonal driver for S^2
      void compute_diagonal_spin_block(DSubSpace& subspace, std::map<std::pair<int, int>, double>& spinmap);

      // Off-diagonal stuff
      std::shared_ptr<Matrix> couple_blocks(DSubSpace& AB, DSubSpace& ApBp); // Off-diagonal driver for H

      std::shared_ptr<Matrix> compute_inter_2e(DSubSpace& AB, DSubSpace& ApBp);
      std::shared_ptr<Matrix> compute_aET(DSubSpace& AB, DSubSpace& ApBp);
      std::shared_ptr<Matrix> compute_bET(DSubSpace& AB, DSubSpace& ApBp);
      std::shared_ptr<Matrix> compute_abFlip(DSubSpace& AB, DSubSpace& ApBp);
      std::shared_ptr<Matrix> compute_abET(DSubSpace& AB, DSubSpace& ApBp);
      std::shared_ptr<Matrix> compute_aaET(DSubSpace& AB, DSubSpace& ApBp);
      std::shared_ptr<Matrix> compute_bbET(DSubSpace& AB, DSubSpace& ApBp);
};

// Locks to make sure the following files are not included on their own
#define MEH_HEADERS
#include <src/meh/meh_compute.hpp>
#include <src/meh/meh_compute_diagonal.hpp>
#include <src/meh/meh_compute_offdiagonal.hpp>
#include <src/meh/meh_gamma_coupling.hpp>
#include <src/meh/meh_init.hpp>
#include <src/meh/meh_spin_coupling.hpp>
#undef MEH_HEADERS

}

#endif
