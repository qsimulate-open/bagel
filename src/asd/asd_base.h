//
// BAGEL - Parallel electron correlation program.
// Filename: asd/asd_base.h
// Copyright (C) 2014 Shane Parker
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

    double thresh_;
    double print_thresh_;

    bool fixed_ci_;

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

    //RDM debug functions
    void debug_rdm(std::shared_ptr<RDM<1>>&, std::shared_ptr<RDM<2>>&, const int istate, const bool mute) const;
    void debug_energy(std::shared_ptr<RDM<1>>&, std::shared_ptr<RDM<2>>&, const int istate, const bool mute) const;

  public:
    ASD_base(const std::shared_ptr<const PTree> input, std::shared_ptr<const Dimer> dimer);

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

    void update_dimer(std::shared_ptr<const Dimer> dimer);

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



#if 0
    // Dimer reduced density matrices
    std::shared_ptr<RDM<1>> onerdm_; // First-order RDM
    std::shared_ptr<RDM<2>> twordm_; // Second-order RDM
    std::shared_ptr<RDM<2>> approx2rdm_; // Second-order RDM

    std::shared_ptr<RDM<3>> threerdm_; // Third-order RDM
    std::map<std::string,std::shared_ptr<Matrix>> rdm3_;
    std::shared_ptr<RDM<4>> fourrdm_; // Fourth-order RDM
    std::map<std::string,std::shared_ptr<Matrix>> fourrdmparts_;

    //3RDM
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> couple_blocks_RDM34(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_diag_RDM34(const std::array<MonomerKey,4>&, const bool) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abFlip_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bbET_RDM34(const std::array<MonomerKey,4>&) const;

    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aabET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abbET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaaET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bbbET_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aETFlip_RDM34(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bETFlip_RDM34(const std::array<MonomerKey,4>&) const;

    void initialize_3RDM();
    void initialize_4RDM();
    double element_3RDM(const int a, const int i, const int b, const int j,const int c, const int k) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_diagonal_block_RDM34(const DimerSubspace_base& AB) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_diagonal_block_RDM4(const DimerSubspace_base& AB) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> couple_blocks_4RDM(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_diag_4RDM(const std::array<MonomerKey,4>&, const bool) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abFlip_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bbET_4RDM(const std::array<MonomerKey,4>&) const;

    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aabET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abbET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaaET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bbbET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aETFlip_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bETFlip_4RDM(const std::array<MonomerKey,4>&) const;

    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaaaET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaabET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aabbET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_abbbET_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bbbbET_4RDM(const std::array<MonomerKey,4>&) const;

    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_aaETFlip_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_bbETFlip_4RDM(const std::array<MonomerKey,4>&) const;
    std::tuple<std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> compute_doubleFlip_4RDM(const std::array<MonomerKey,4>&) const;

    void symmetrize_RDM34() const;
    void symmetrize_RDM4();


  template<int a,  int i,  int b,  int j,  int c,  int k,  int d,  int l, // according to the sorted
           int Xa, int Xi, int Xb, int Xj, int Xc, int Xk, int Xd, int Xl, // 0 for A, 1 for B
           bool /*transpose*/trans >
  void fill_RDM(std::shared_ptr<Matrix> rdmt,
                const int n0, const int n1, const int n2, const int n3, const int n4, const int n5, const int n6, const int n7) { // according to the unsorted

    const int nactA = dimer_->embedded_refs().first->nact();
    const int nactB = dimer_->embedded_refs().second->nact();
    const int nactT = nactA + nactB;

    //                  A       B
    int la = Xa == 0 ? 0     : nactA;
    int li = Xi == 0 ? 0     : nactA;
    int lb = Xb == 0 ? 0     : nactA;
    int lj = Xj == 0 ? 0     : nactA;
    int lc = Xc == 0 ? 0     : nactA;
    int lk = Xk == 0 ? 0     : nactA;
    int ld = Xd == 0 ? 0     : nactA;
    int ll = Xl == 0 ? 0     : nactA;

    int ha = Xa == 0 ? nactA : nactT;
    int hi = Xi == 0 ? nactA : nactT;
    int hb = Xb == 0 ? nactA : nactT;
    int hj = Xj == 0 ? nactA : nactT;
    int hc = Xc == 0 ? nactA : nactT;
    int hk = Xk == 0 ? nactA : nactT;
    int hd = Xd == 0 ? nactA : nactT;
    int hl = Xl == 0 ? nactA : nactT;

    auto mat = std::make_shared<Matrix>(n0*n1*n2*n3*n4*n5*n6*n7,1); //empty

    if (trans) { //transposed
                        //sorted                                                unsorted dim/s
      sort_indices<i,a,j,b,k,c,l,d, 0,1, 1,1>(rdmt->data(), mat->data(), n0, n1, n2, n3, n4, n5, n6, n7);
      auto low = {li,la,lj,lb,lk,lc,ll,ld};
      auto up  = {hi,ha,hj,hb,hk,hc,hl,hd};
      auto outv = btas::make_rwview(fourrdm_->range().slice(low,up), fourrdm_->storage());
      copy(mat->begin(), mat->end(), outv.begin());
    } else { //normal
                        //sorted                                                unsorted dim/s
      sort_indices<a,i,b,j,c,k,d,l, 0,1, 1,1>(rdmt->data(), mat->data(), n0, n1, n2, n3, n4, n5, n6, n7);
      auto low = {la,li,lb,lj,lc,lk,ld,ll};
      auto up  = {ha,hi,hb,hj,hc,hk,hd,hl};
      auto outv = btas::make_rwview(fourrdm_->range().slice(low,up), fourrdm_->storage());
      copy(mat->begin(), mat->end(), outv.begin());
    }
  }

#include <src/util/prim_op.h>

#include <src/wfn/reference.h>

#endif
