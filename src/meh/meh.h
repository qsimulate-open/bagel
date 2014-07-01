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
#include <src/meh/meh_base.h>
#include <src/smith/prim_op.h>

namespace bagel {

/// Wrapper for monomer CI wavefunctions that includes extra helpful information
template <class VecType>
struct MonomerSubspace_base {
  const int S_;
  const int ms_;
  const int charge_;
  std::shared_ptr<const VecType> monomerci_;

  MonomerSubspace_base(const int S, const int ms, const int charge, std::shared_ptr<const VecType> monomerci) :
    S_(S), ms_(ms), charge_(charge), monomerci_(monomerci) {}
};

/// Contains all of the information for a product of two monomer spaces
template <class VecType>
class DimerSubspace_base {
  protected:
    const int offset_;
    const int nstatesA_;
    const int nstatesB_;
    const std::string stringA_;
    const std::string stringB_;

    std::pair<MonomerSubspace_base<VecType>, MonomerSubspace_base<VecType>> ci_;
    std::pair<std::shared_ptr<const CSymMatrix>, std::shared_ptr<const CSymMatrix>> sigma_;

  public:
    DimerSubspace_base(int& _offset, const SpaceKey Akey, const SpaceKey Bkey, std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> _ci) :
      offset_(_offset), nstatesA_(_ci.first->ij()), nstatesB_(_ci.second->ij()), stringA_(Akey.to_string()), stringB_(Bkey.to_string()),
       ci_(std::make_pair(MonomerSubspace_base<VecType>(Akey.S, Akey.m_s, Akey.q, _ci.first),
                     MonomerSubspace_base<VecType>(Bkey.S, Bkey.m_s, Bkey.q, _ci.second)))
    { _offset += dimerstates(); }

    const int offset() const { return offset_; }
    const int dimerstates() const { return nstatesA_ * nstatesB_; }
    const int dimerindex(const int iA, const int iB) const { return (iA + iB*nstatesA_); }
    const std::string string(const int i, const int j) const {
      std::string out = stringA_ + lexical_cast<std::string>(i) + std::string(" ") + stringB_ + lexical_cast<std::string>(j);
      return out;
    }

    template <int unit> const int nstates() const { return ( unit == 0 ? nstatesA_ : nstatesB_ ); }
    template <int unit> std::shared_ptr<const VecType> ci() const { return ( unit == 0 ? ci_.first.monomerci_ : ci_.second.monomerci_ ); }
    template <int unit> std::shared_ptr<const CSymMatrix> sigma() const { return ( unit == 0 ? sigma_.first : sigma_.second ); }

    std::pair<int,int> S() const { return std::make_pair(ci_.first.S_, ci_.second.S_); }
    std::pair<int,int> ms() const { return std::make_pair(ci_.first.ms_, ci_.second.ms_); }
    std::pair<int,int> charge() const { return std::make_pair(ci_.first.charge_, ci_.second.charge_); }

    template <int unit> void set_sigma(std::shared_ptr<const CSymMatrix> s) { (unit == 0 ? sigma_.first : sigma_.second) = s; }

};


// forward declaration
template <class VecType>
class MultiExcitonHamiltonian;

namespace asd {
// implementation class. true = Hamiltonian, false = RDM.
template <bool _N>
struct ASD_impl {
  template<typename T>
  using DSb = DimerSubspace_base<T>;
  using return_type = typename std::conditional<_N, Matrix, RDM<2>>::type;
  template <class VecType>
  static std::shared_ptr<return_type> compute_diagonal_block(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& subspace) { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_inter_2e(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp) { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_aET(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp)  { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_bET(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp)  { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_abFlip(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp) { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_abET(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp) { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_aaET(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp) { assert(false); return nullptr; }
  template <class VecType>
  static std::shared_ptr<return_type> compute_bbET(MultiExcitonHamiltonian<VecType>*, DSb<VecType>& AB, DSb<VecType>& ApBp) { assert(false); return nullptr; }
};
}

/// Template for MEH (to be renamed ASD)
template <class VecType>
class MultiExcitonHamiltonian : public MEH_base {
  protected: using DSubSpace = DimerSubspace_base<VecType>;
  protected: using DCISpace = DimerCISpace_base<VecType>;
  protected: using CiType = typename VecType::Ci;

  template <bool>
  friend class asd::ASD_impl;

  protected:
    std::shared_ptr<DCISpace> cispace_;
    std::shared_ptr<GammaForest<VecType, 2>> gammaforest_;
    std::vector<DSubSpace> subspaces_;

    void modelize();

  public:
    MultiExcitonHamiltonian(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DCISpace> cispace);

    void compute() override;
    void compute_rdm() const override;

    void print_hamiltonian(const std::string title = "MultiExciton Hamiltonian", const int nstates = 10) const;
    void print_states(const Matrix& cc, const std::vector<double>& energies, const double thresh = 0.01, const std::string title = "Adiabats") const;
    void print_property(const std::string label, std::shared_ptr<const Matrix>, const int size = 10) const ;
    void print(const double thresh = 0.01) const;

  private:
    const Coupling coupling_type(const DSubSpace& AB, const DSubSpace& ApBp) const;

    void generate_initial_guess(std::shared_ptr<Matrix> cc, std::vector<DSubSpace>& subspace, const int nstates);
    std::shared_ptr<Matrix> apply_hamiltonian(const Matrix& o, std::vector<DSubSpace>& subspaces);
    std::vector<double> diagonalize(std::shared_ptr<Matrix>& cc, std::vector<DSubSpace>& subspace, const bool mute = false);

    std::shared_ptr<Matrix> compute_1e_prop(std::shared_ptr<const Matrix> hAA, std::shared_ptr<const Matrix> hBB, std::shared_ptr<const Matrix> hAB, const double core) const;
    std::shared_ptr<Matrix> compute_offdiagonal_1e(const DSubSpace& AB, const DSubSpace& ApBp, std::shared_ptr<const Matrix> hAB) const;
    std::shared_ptr<Matrix> compute_diagonal_1e(const DSubSpace& subspace, const double* hAA, const double* hBB, const double diag) const;

    // Diagonal block stuff
    void compute_pure_terms(DSubSpace& subspace, std::shared_ptr<const DimerJop> jop);
    std::shared_ptr<Matrix> compute_intra(const DSubSpace& subspace, std::shared_ptr<const DimerJop> jop, const double diag);

    virtual std::shared_ptr<VecType> form_sigma(std::shared_ptr<const VecType> ccvec, std::shared_ptr<const MOFile> jop) const = 0;
    virtual std::shared_ptr<VecType> form_sigma_1e(std::shared_ptr<const VecType> ccvec, const double* modata) const = 0;

    // Gamma Tree building
    void gamma_couple_blocks(DSubSpace& AB, DSubSpace& ApBp);
    void spin_couple_blocks(DSubSpace& AB, DSubSpace& ApBp, std::map<std::pair<int, int>, double>& spinmap); // Off-diagonal driver for S^2
    void compute_diagonal_spin_block(DSubSpace& subspace, std::map<std::pair<int, int>, double>& spinmap);

    // Off-diagonal stuff
    template <bool _N, typename return_type = typename std::conditional<_N, Matrix, RDM<2>>::type>
    std::shared_ptr<return_type> couple_blocks(DSubSpace& AB, DSubSpace& ApBp); // Off-diagonal driver for H

    template <bool _N> std::shared_ptr<Matrix> compute_diagonal_block(DSubSpace& subspace)      { return asd::ASD_impl<_N>::compute_diagonal_block(this, subspace); }
    template <bool _N> std::shared_ptr<Matrix> compute_inter_2e(DSubSpace& AB, DSubSpace& ApBp) { return asd::ASD_impl<_N>::compute_inter_2e(this, AB, ApBp); }
    template <bool _N> std::shared_ptr<Matrix> compute_aET(DSubSpace& AB, DSubSpace& ApBp)      { return asd::ASD_impl<_N>::compute_aET(this, AB, ApBp); }
    template <bool _N> std::shared_ptr<Matrix> compute_bET(DSubSpace& AB, DSubSpace& ApBp)      { return asd::ASD_impl<_N>::compute_bET(this, AB, ApBp); }
    template <bool _N> std::shared_ptr<Matrix> compute_abFlip(DSubSpace& AB, DSubSpace& ApBp)   { return asd::ASD_impl<_N>::compute_abFlip(this, AB, ApBp); }
    template <bool _N> std::shared_ptr<Matrix> compute_abET(DSubSpace& AB, DSubSpace& ApBp)     { return asd::ASD_impl<_N>::compute_abET(this, AB, ApBp); }
    template <bool _N> std::shared_ptr<Matrix> compute_aaET(DSubSpace& AB, DSubSpace& ApBp)     { return asd::ASD_impl<_N>::compute_aaET(this, AB, ApBp); }
    template <bool _N> std::shared_ptr<Matrix> compute_bbET(DSubSpace& AB, DSubSpace& ApBp)     { return asd::ASD_impl<_N>::compute_bbET(this, AB, ApBp); }
};

// Locks to make sure the following files are not included on their own
#define MEH_HEADERS
#include <src/meh/meh_compute.hpp>
#include <src/meh/meh_compute_diagonal.hpp>
#include <src/meh/meh_compute_offdiagonal.hpp>
#include <src/meh/meh_gamma_coupling.hpp>
#include <src/meh/meh_init.hpp>
#include <src/meh/meh_modelize.hpp>
#include <src/meh/meh_diagonalize.hpp>
#include <src/meh/meh_spin_coupling.hpp>
#include <src/meh/meh_compute_rdm.hpp>
#undef MEH_HEADERS

}

#endif
