//
// BAGEL - Parallel electron correlation program.
// Filename: zfci.h
// Copyright (C) 2013 Toru Shiozaki 
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu> 
// Modified by: Michael Caldwell  <caldwell@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#ifndef __BAGEL_ZFCI_ZFCI_H
#define __BAGEL_ZFCI_ZFCI_H

#include <src/fci/space.h>
#include <src/fci/dvec.h>
#include <src/zfci/zmofile.h>
#include <src/zfci/relmofile.h>
#include <src/wfn/ciwfn.h>
#include <src/math/zmatrix.h>

namespace bagel {

class ZFCI : public Method {

  protected:
    // max #iteration
    int max_iter_;
    // threshold for variants
    double thresh_;
    double print_thresh_;

    // numbers of electrons
    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;

    // number of states
    int nstate_;

#if 0
    // properties to be calculated
    std::vector<std::shared_ptr<CIProperties>> properties_;
#endif

    // total energy
    std::vector<double> energy_;

    // CI vector at convergence
    std::shared_ptr<ZDvec> cc_;
    // RDMs; should be resized in constructors
    std::vector<std::shared_ptr<RDM<1>>> rdm1_;
    std::vector<std::shared_ptr<RDM<2>>> rdm2_;
    // state averaged RDM
    std::vector<double> weight_;
    std::shared_ptr<RDM<1>> rdm1_av_;
    std::shared_ptr<RDM<2>> rdm2_av_;
    // MO integrals
    std::shared_ptr<ZMOFile_Base> jop_;


    // Determinant space
    std::shared_ptr<const Space_base> space_;

    // denominator
    std::shared_ptr<Civec> denom_;

    // some init functions
    void common_init(); // may end up unnecessary
    // obtain determinants for guess generation
    void generate_guess(const int nspin, const int nstate, std::shared_ptr<ZDvec>);
    // generate spin-adapted guess configurations
    virtual std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet);

    /* Virtual functions -- these MUST be defined in the derived class*/
    // denominator
    virtual void const_denom() = 0;
#if 0
    // functions related to natural orbitals
    void update_rdms(const std::shared_ptr<ZMatrix>& coeff);

    // internal function for RDM1 and RDM2 computations
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_last_step(std::shared_ptr<const ZDvec>, std::shared_ptr<const ZDvec>, std::shared_ptr<const ZCivec>) const;
#endif

    // print functions
    void print_header() const;

  public:
    // this constructor is ugly... to be fixed some day...
    ZFCI(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>,
         const int ncore = -1, const int nocc = -1, const int nstate = -1);

    virtual void compute() override;

    virtual void update(std::shared_ptr<const Coeff> ) = 0;

    // returns members
    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    double core_energy() const { return jop_->core_energy(); }

    virtual int nij() const { return norb_*norb_; }

    double weight(const int i) const { return weight_[i]; }

    // virtual application of Hamiltonian
    virtual std::shared_ptr<ZDvec> form_sigma(std::shared_ptr<const ZDvec> c, std::shared_ptr<const ZMOFile_Base> jop, const std::vector<int>& conv) const = 0;


#if 0
    // rdms
    void compute_rdm12(); // compute all states at once + averaged rdm
    void compute_rdm12(const int istate);
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> compute_rdm34(const int istate) const;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_from_civec(std::shared_ptr<const ZCivec>, std::shared_ptr<const ZCivec>) const;
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_av_from_dvec(std::shared_ptr<const ZDvec>, std::shared_ptr<const ZDvec>,
                                 std::shared_ptr<const Determinants> o = std::shared_ptr<const Determinants>()) const;

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
    // move to natural orbitals
    std::pair<std::shared_ptr<ZMatrix>,std::vector<std::complex<double>>> natorb_convert();
#endif

    std::shared_ptr<const Space_base> space() const { return space_; }

    // returns integral files
    std::shared_ptr<const ZMOFile_Base> jop() const { return jop_; }

    // returns a denominator
    //returns non-complex because denom_ is currently all real because of real mo2e
    std::shared_ptr<const Civec> denom() const { return denom_; }

    // returns total energy
    std::vector<double> energy() const { return energy_; }
    double energy(const int i) const { return energy_.at(i); }

    // returns CI vectors
    std::shared_ptr<ZDvec> civectors() const { return cc_; }

    // These are needed for the RDM stuff, apparently
    void sigma_2a1(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;
    void sigma_2a2(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;

    std::shared_ptr<const CIWfn> conv_to_ciwfn();
    std::shared_ptr<const Reference> conv_to_ref() const override { return std::shared_ptr<const Reference>(); }
};

}

#endif

