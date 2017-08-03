//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/rasci.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#ifndef __BAGEL_RAS_RASCI_H
#define __BAGEL_RAS_RASCI_H

#define NORDMS

#include <src/ci/fci/mofile.h>
#include <src/ci/ras/civector.h>
#include <src/wfn/method.h>
#include <src/wfn/reference.h>

namespace bagel {

class RASCI : public Method {
  protected:
    // max #iteration
    int max_iter_;
    int davidson_subspace_;
    int nguess_;

    // threshold for variants
    double thresh_;
    double print_thresh_;

    // algorithmic options
    int batchsize_;

    // numbers of electrons
    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;
    std::array<int,3> ras_;
    int max_holes_;
    int max_particles_;

    // number of states
    int nstate_;

    // total energy
    std::vector<double> energy_;

    // CI vector at convergence
    std::shared_ptr<RASDvec> cc_;

#ifndef NORDMS
    // RDMs; should be resized in constructors
    std::vector<std::shared_ptr<RDM<1>>> rdm1_;
    std::vector<std::shared_ptr<RDM<2>>> rdm2_;
    // state averaged RDM
    std::vector<double> weight_;
    std::shared_ptr<RDM<1>> rdm1_av_;
    std::shared_ptr<RDM<2>> rdm2_av_;
#endif

    // MO integrals
    std::shared_ptr<MOFile> jop_;

    // Determinant space
    std::shared_ptr<const RASDeterminants> det_;

    // denominator
    std::shared_ptr<RASCivec> denom_;

    // some init functions
    void common_init(); // may end up unnecessary
    void create_Jiiii();
    // obtain determinants for guess generation
    void generate_guess(const int nspin, const int nstate, std::shared_ptr<RASDvec>& out);
    void model_guess(std::shared_ptr<RASDvec>& out);
    // generate spin-adapted guess configurations
    std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet);

    // denominator
    void const_denom();

    // functions related to natural orbitals
    void update_rdms(const std::shared_ptr<Matrix>& coeff);

#ifndef NORDMS
    // internal function for RDM1 and RDM2 computations
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_last_step(std::shared_ptr<RASDvec>, std::shared_ptr<const RASDvec>, std::shared_ptr<const RASCivec>) const;
#endif

    // print functions
    void print_header() const;

  public:
    // this constructor is ugly... to be fixed some day...
    RASCI(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;

    void update(std::shared_ptr<const Coeff>);

    // returns members
    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    template <int space> int ras() const { return std::get<space>(ras_); }
    double core_energy() const { return jop_->core_energy(); }

#ifndef NORDMS
    // rdms
    void compute_rdm12(); // compute all states at once + averaged rdm
    void compute_rdm12(const int istate);
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> compute_rdm34(const int istate) const;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_from_civec(std::shared_ptr<const RASCivec>, std::shared_ptr<const RASCivec>) const;
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_av_from_dvec(std::shared_ptr<const RASDvec>, std::shared_ptr<const RASDvec>, std::shared_ptr<const Determinants> o = nullptr) const;

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

    // rdm ci derivatives
    std::shared_ptr<RASDvec> rdm1deriv() const;
    std::shared_ptr<RASDvec> rdm2deriv() const;

    // move to natural orbitals
    std::pair<std::shared_ptr<Matrix>, std::vector<double>> natorb_convert();
#endif

    const std::shared_ptr<const Geometry> geom() const { return geom_; }

    std::shared_ptr<const RASDeterminants> det() const { return det_; }

    // returns integral files
    std::shared_ptr<const MOFile> jop() const { return jop_; }

    // returns a denominator
    std::shared_ptr<const RASCivec> denom() const { return denom_; }

    // returns total energy
    std::vector<double> energy() const { return energy_; }
    double energy(const int i) const { return energy_.at(i); }

    // returns CI vectors
    std::shared_ptr<RASDvec> civectors() const { return cc_; }

    std::shared_ptr<const Reference> conv_to_ref() const override { return nullptr; }
};

}

#endif

