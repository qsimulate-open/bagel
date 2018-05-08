//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_base.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __ASD_ASD_BASE_H
#define __ASD_ASD_BASE_H

#include <src/asd/dimer/dimer.h>
#include <src/asd/dimer/dimer_jop.h>
#include <src/asd/asd_spin.h>
#include <src/asd/gamma_tensor.h>
#include <src/asd/coupling.h>
#include <src/asd/state_tensor.h>

namespace bagel {

/// Specifies a single block of a model Hamiltonian
struct ModelBlock {
  std::pair<int, int> S_;
  std::pair<int, int> charge_;
  std::pair<int, int> M_;
  int nstates_;

  ModelBlock(std::pair<int, int> S, std::pair<int, int> q, std::pair<int, int> M, const int nstates) :
    S_(S), charge_(q), M_(M), nstates_(nstates) {}
};

/// Base class for ASD. Contains everything that doesn't need to be templated.
class ASD_base {
  protected:
    std::shared_ptr<const Dimer> dimer_;

    std::shared_ptr<DimerJop> jop_;

    std::shared_ptr<Matrix> hamiltonian_; ///< if stored_ is true, Hamiltonian is stored here
    std::shared_ptr<Matrix> adiabats_; ///< Eigenvectors of adiabatic states
    std::vector<std::pair<std::string, std::shared_ptr<Matrix>>> properties_;

    std::unique_ptr<double[]> denom_; ///< denominator used for Davidson

    int max_spin_;
    std::shared_ptr<ASDSpin> spin_; ///< Sparsely stored spin matrix

    std::vector<double> energies_; ///< Adiabatic energies

    // RDMs; should be resized in constructors
    std::vector<std::shared_ptr<RDM<1>>> rdm1_;
    std::vector<std::shared_ptr<RDM<2>>> rdm2_;
    // state averaged RDM
    std::vector<double> weight_;
    std::shared_ptr<RDM<1>> rdm1_av_;
    std::shared_ptr<RDM<2>> rdm2_av_;

    // Total system quantities
    int dimerstates_; ///< Total size of dimer Hamiltonian. Counted up during initialization

    // Options
    int nstates_;
    int nspin_;
    int charge_;
    int max_iter_;
    int nguess_;
    int davidson_subspace_;

    bool store_matrix_;
    bool dipoles_;
    bool print_info_;

    double thresh_;
    double print_thresh_;

    // Orbital optimization related
    bool compute_rdm_;
    bool fix_ci_;

    virtual std::vector<DimerSubspace_base> subspaces_base() const = 0;

    // Gamma Tensor
    std::array<std::shared_ptr<const GammaTensor>,2> gammatensor_;

    std::shared_ptr<StateTensor> statetensor_;

    std::vector<std::vector<ModelBlock>> models_to_form_; ///< Contains specifications to construct model spaces
    std::vector<std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>>> models_; ///< models that have been built

    std::shared_ptr<Matrix> apply_hamiltonian(const Matrix& o, const std::vector<DimerSubspace_base>& subspaces);
    std::vector<double> diagonalize(std::shared_ptr<Matrix>& cc, const std::vector<DimerSubspace_base>& subspace, const bool mute = false);

    // Off-diagonal stuff
    std::shared_ptr<Matrix> couple_blocks(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const; // Off-diagonal driver for H

    std::shared_ptr<Matrix> compute_offdiagonal_1e(const std::array<MonomerKey,4>&, std::shared_ptr<const Matrix> h) const;
    std::shared_ptr<Matrix> compute_inter_2e(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_aET(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_bET(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_abFlip(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_abET(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_aaET(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_bbET(const std::array<MonomerKey,4>&) const;
    std::shared_ptr<Matrix> compute_diagonal_block(const DimerSubspace_base& subspace) const;

    void generate_initial_guess(std::shared_ptr<Matrix> cc, const std::vector<DimerSubspace_base>& subspace, const int nstates);
    std::shared_ptr<Matrix> compute_intra(const DimerSubspace_base& subspace, std::shared_ptr<const DimerJop> jop, const double diag) const;

    void modelize();

    void print_hamiltonian(const std::string title = "MultiExciton Hamiltonian", const int nstates = 10) const;
    void print_states(const Matrix& cc, const std::vector<double>& energies, const double thresh = 0.01, const std::string title = "Adiabats") const;
    void print_property(const std::string label, std::shared_ptr<const Matrix>, const int size = 10) const ;
    void print(const double thresh = 0.01) const;

    void compute_rdm12_dimer();

    //RDM information print
    void print_rdm_info(std::shared_ptr<RDM<1>>&, std::shared_ptr<RDM<2>>&, const int istate) const;
    void print_energy_info(std::shared_ptr<RDM<1>>&, std::shared_ptr<RDM<2>>&, const int istate) const;

  public:
    ASD_base(const std::shared_ptr<const PTree> input, std::shared_ptr<const Dimer> dimer, bool rdm = false);

    virtual void compute() = 0;

    std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> model(const int i) { return models_[i]; }

    std::vector<double> energy() const { return energies_; }
    double energy(const int i) const { return energies_.at(i); }

    std::vector<std::shared_ptr<RDM<1>>> rdm1() { return rdm1_; }
    std::vector<std::shared_ptr<RDM<2>>> rdm2() { return rdm2_; }
    std::shared_ptr<RDM<1>> rdm1(const int i) { return rdm1_.at(i); }
    std::shared_ptr<RDM<2>> rdm2(const int i) { return rdm2_.at(i); }
    std::shared_ptr<const RDM<1>> rdm1(const int i) const { return rdm1_.at(i); }
    std::shared_ptr<const RDM<2>> rdm2(const int i) const { return rdm2_.at(i); }
    std::shared_ptr<RDM<1>> rdm1_av() { return rdm1_av_; }
    std::shared_ptr<RDM<2>> rdm2_av() { return rdm2_av_; }
    std::shared_ptr<const RDM<1>> rdm1_av() const { return rdm1_av_; }
    std::shared_ptr<const RDM<2>> rdm2_av() const { return rdm2_av_; }

    void update_dimer_and_fix_ci(std::shared_ptr<const Dimer> dimer);

  private:
    // RDM
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_diagonal_block(const DimerSubspace_base& subspace, const int istate) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> couple_blocks(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp, const int istate) const;

    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_inter_2e(const std::array<MonomerKey,4>&, const int istate, const bool subdia) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_aET(const std::array<MonomerKey,4>&, const int istate) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_bET(const std::array<MonomerKey,4>&, const int istate) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_abFlip(const std::array<MonomerKey,4>&, const int istate) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_abET(const std::array<MonomerKey,4>&, const int istate) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_aaET(const std::array<MonomerKey,4>&, const int istate) const;
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_bbET(const std::array<MonomerKey,4>&, const int istate) const;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>> compute_rdm12_dimer(const int i) const;
    void symmetrize_rdm12(std::shared_ptr<RDM<1>>&, std::shared_ptr<RDM<2>>&) const;
};

}

#endif
