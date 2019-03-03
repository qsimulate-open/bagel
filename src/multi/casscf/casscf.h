//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casscf.h
// Copyright (C) 2011 Toru Shiozaki
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

// This is a base class for various CASSCF solvers.
// The assumpation is made that the CI and orbital optimizations are done independently.
// This should be a good strategy in large systems.
//

#ifndef __BAGEL_CASSCF_CASSCF_H
#define __BAGEL_CASSCF_CASSCF_H

#include <src/wfn/reference.h>
#include <src/util/muffle.h>
#include <src/ci/fci/distfci.h>
#include <src/ci/fci/knowles.h>
#include <src/ci/fci/harrison.h>
#include <src/multi/casscf/rotfile.h>

namespace bagel {

enum FCIAlgorithmType { knowles, harrison, dist };

class FCI_algorithms {
  protected:
    FCIAlgorithmType type_;

  public:
    FCI_algorithms() : type_(FCIAlgorithmType::knowles) { }
    FCI_algorithms(std::string input_algorithm) {
      if (input_algorithm == "knowles" || input_algorithm == "kh" || input_algorithm == "handy") {
        type_ = FCIAlgorithmType::knowles;
      } else if (input_algorithm == "harrison" || input_algorithm == "zarrabian" || input_algorithm == "hz") {
        type_ = FCIAlgorithmType::harrison;
#ifdef HAVE_MPI_H
      } else if (input_algorithm == "parallel" || input_algorithm == "dist") {
        type_ = FCIAlgorithmType::dist;
#endif
      } else {
        throw std::runtime_error("Unknown FCI algorithm specified. " + input_algorithm);
      }
    }

    bool is_knowles() const { return type_ == FCIAlgorithmType::knowles; }
    bool is_harrison() const { return type_ == FCIAlgorithmType::harrison; }
    bool is_dist() const { return type_ == FCIAlgorithmType::dist; }
};


class CASSCF : public Method, public std::enable_shared_from_this<CASSCF> {

  protected:
    // some internal information
    int nocc_; // sum of nact_ + nclosed_
    int nclosed_;
    int nact_;
    int nvirt_;
    int nmo_;
    int nstate_;
    int max_iter_;
    int max_micro_iter_;
    double thresh_;
    double thresh_micro_;
    double thresh_overlap_;
    bool conv_ignore_;
    bool restart_cas_;
    bool natocc_;
    bool canonical_;

    std::shared_ptr<const Coeff> coeff_;

    // RDMs are given externally (e.g., FCIQMC)
    std::string external_rdm_;
    std::shared_ptr<FCI_algorithms> fci_algorithm_;
    std::shared_ptr<FCI_base> fci_;
    void print_header() const;
    void print_iteration(const int iter, const std::vector<double>& energy, const double error, const double time) const;
    void common_init();

    void mute_stdcout();
    void resume_stdcout();

    const std::shared_ptr<const Matrix> hcore_;

    std::tuple<std::shared_ptr<const Coeff>,VectorB,VectorB> semi_canonical_orb() const;
    std::shared_ptr<const Matrix> spin_density() const;

    // orbital eigenvalues and occupations
    VectorB eig_;
    VectorB occup_;

    // energy
    std::vector<double> energy_;
    double rms_grad_;

    // properties
    bool do_hyperfine_;

    // mask some of the output
    mutable std::shared_ptr<Muffle> muffle_;

  public:
    CASSCF(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> = nullptr);
    virtual ~CASSCF();

    virtual void compute() override = 0;

    std::shared_ptr<const Reference> ref() const { return ref_; }
    virtual std::shared_ptr<const Reference> conv_to_ref() const override;

    std::shared_ptr<FCI_base> fci() { return fci_; }

    // functions to retrieve protected members
    int nocc() const { return nocc_; }
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nvirt() const { return nvirt_; }
    int nmo() const { return nmo_; }
    int nstate() const { return nstate_; }
    int max_iter() const { return max_iter_; }
    int max_micro_iter() const { return max_micro_iter_; }
    double thresh() const { return thresh_; }
    double thresh_micro() const { return thresh_micro_; }
    double thresh_overlap() const { return thresh_overlap_; }

    double energy(const int i) const { return energy_[i]; }
    double energy_av() const { return blas::average(energy_); }
    const std::vector<double>& energy() const { return energy_; }
    double rms_grad() const { return rms_grad_; }

    std::shared_ptr<Matrix> compute_active_fock(const MatView acoeff, std::shared_ptr<const RDM<1>> rdm1) const;

    std::shared_ptr<Matrix> ao_rdm1(std::shared_ptr<const RDM<1>> rdm1, const bool inactive_only = false) const;
    std::shared_ptr<const Matrix> hcore() const { return hcore_; }

    std::shared_ptr<const Coeff> coeff() const { return coeff_; }
};

}

#endif
