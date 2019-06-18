//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_base.h
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_CI_FCI_FCI_BASE_H
#define __SRC_CI_FCI_FCI_BASE_H

#include <src/ci/fci/mofile.h>
#include <src/wfn/ciwfn.h>
#include <src/wfn/method.h>
#include <src/ci/fci/dvec.h>
#include <src/ci/fci/distcivec.h>
#include <src/util/math/davidson.h>

namespace bagel {

class FCI_base : public Method {
  protected:

    // Options
    int max_iter_;
    int davidson_subspace_;
    int nguess_;
    double thresh_;
    double print_thresh_;

    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;

    int nstate_;

    // extra
    std::shared_ptr<const Determinants> det_;

    // results
    std::vector<double> energy_;
    std::shared_ptr<MOFile> jop_;
    std::shared_ptr<const Matrix> coeff_;

    // RDMs; should be resized in constructors
    std::shared_ptr<VecRDM<1>> rdm1_;
    std::shared_ptr<VecRDM<2>> rdm2_;
    // state averaged RDM
    std::vector<double> weight_;
    std::shared_ptr<RDM<1>> rdm1_av_;
    std::shared_ptr<RDM<2>> rdm2_av_;

    virtual void print_header() const = 0;

    // const_denom function here only makes a denom for local data of DistCivec.
    virtual void const_denom() = 0;

    // restart
    bool restart_;
    bool restarted_;

    // integral reuse
    bool store_half_ints_;

    // functions related to natural orbitals
    void update_rdms(std::shared_ptr<const Matrix> coeff);

  public:
    // this constructor is ugly... to be fixed some day...
    FCI_base(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r,
             const int ncore = -1, const int norb = -1, const int nstate = -1, const bool store = false)
      : Method(idat, g, r), ncore_(ncore), norb_(norb), nstate_(nstate), restarted_(false), store_half_ints_(store) {
    }

    FCI_base() { }
    virtual ~FCI_base() { }

    // FCI compute function
    virtual void compute() override = 0;

    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    double core_energy() const { return jop_->core_energy(); }
    double weight(const int i) const { return weight_.at(i); }

    virtual void update(std::shared_ptr<const Matrix>) = 0;

    std::shared_ptr<const Determinants> det() const { return det_; }
    std::shared_ptr<const MOFile> jop() const { return jop_; }
    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
    virtual std::shared_ptr<const Civec> denom() const = 0;
    virtual std::shared_ptr<const Dvec> civectors() const = 0;

    std::vector<double> energy() const { return energy_; }
    double energy(const int i) const { return energy_.at(i); }

    // compute all states at once + averaged rdm
    virtual void compute_rdm12() = 0;
    virtual void compute_rdm12(const int ist, const int jst) = 0;
    // compute 3 and 4 RDMs
    virtual std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> rdm34(const int ist, const int jst) const = 0;
    virtual std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<3>>> rdm34f(const int ist, const int jst, std::shared_ptr<const Matrix> fock) const = 0;
    // compute "alpha" 1 and 2 RDMs <ia ja> and <ia ja, k, l>
    virtual std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>> rdm12_alpha(const int ist, const int jst) const = 0;
    // compute "alpha" 3 and 4 RDMs <ia ja, k, l, m n>...
    virtual std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> rdm34_alpha(const int ist, const int jst) const = 0;

    virtual std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_av_from_dvec(std::shared_ptr<const Dvec>, std::shared_ptr<const Dvec>, std::shared_ptr<const Determinants> o = nullptr) const = 0;

    // rdm ci derivatives
    virtual std::shared_ptr<Dvec> rdm1deriv(const int istate) const = 0;
    virtual std::shared_ptr<Dvec> rdm2deriv(const int istate) const = 0;
    virtual std::shared_ptr<Matrix> rdm2fderiv(const int istate, std::shared_ptr<const Matrix> fock, std::shared_ptr<const Matrix> dmat) const = 0;
    virtual std::shared_ptr<Matrix> rdm2deriv_offset(const int istate, const size_t dsize, const size_t offset, std::shared_ptr<const Matrix> dmat, const bool parallel = true) const = 0;
    virtual std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>
      rdm3deriv(const int istate, std::shared_ptr<const Matrix> fock, const size_t offset, const size_t size, std::shared_ptr<const Matrix> dbra_in, std::shared_ptr<const Matrix> fock_ebra_in) const = 0;

    std::shared_ptr<VecRDM<1>> rdm1() { return rdm1_; }
    std::shared_ptr<VecRDM<2>> rdm2() { return rdm2_; }
    std::shared_ptr<RDM<1>> rdm1(const int i) { return rdm1(i, i); }
    std::shared_ptr<RDM<2>> rdm2(const int i) { return rdm2(i, i); }
    std::shared_ptr<RDM<1>> rdm1(const int i, const int j) { return rdm1_->at(i, j); }
    std::shared_ptr<RDM<2>> rdm2(const int i, const int j) { return rdm2_->at(i, j); }
    std::shared_ptr<const RDM<1>> rdm1(const int i) const { return rdm1(i, i); }
    std::shared_ptr<const RDM<2>> rdm2(const int i) const { return rdm2(i, i); }
    std::shared_ptr<const RDM<1>> rdm1(const int i, const int j) const { return rdm1_->at(i, j); }
    std::shared_ptr<const RDM<2>> rdm2(const int i, const int j) const { return rdm2_->at(i, j); }
    std::shared_ptr<RDM<1>> rdm1_av() { return rdm1_av_; }
    std::shared_ptr<RDM<2>> rdm2_av() { return rdm2_av_; }
    std::shared_ptr<const RDM<1>> rdm1_av() const { return rdm1_av_; }
    std::shared_ptr<const RDM<2>> rdm2_av() const { return rdm2_av_; }

    void rotate_rdms(std::shared_ptr<const Matrix> trans);

    virtual std::shared_ptr<const CIWfn> conv_to_ciwfn() const = 0;
    virtual std::shared_ptr<const Reference> conv_to_ref() const = 0;

    virtual void read_external_rdm12_av(const std::string& file) = 0;
};

}

#endif
