//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasscf.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_ZCASSCF_ZCASSCF_H
#define __SRC_ZCASSCF_ZCASSCF_H

#include <src/ci/zfci/zharrison.h>
#include <src/multi/casscf/rotfile.h>
#include <src/util/muffle.h>

namespace bagel {

class ZCASSCF : public Method, public std::enable_shared_from_this<ZCASSCF> {
  protected:
    int nneg_;
    int nocc_;
    int nclosed_;
    int nact_;
    int nvirt_;
    int nvirtnr_;
    int nbasis_;

    int charge_;

    bool gaunt_;
    bool breit_;
    bool natocc_;
    bool canonical_;

    // set if RDMs are given externally (e.g., FCIQMC)
    std::string external_rdm_;

    double thresh_;
    double thresh_micro_;
    double thresh_overlap_;
    bool conv_ignore_;
    bool restart_cas_;

    int nstate_;

    int max_iter_;
    int max_micro_iter_;

    std::shared_ptr<const ZCoeff_Block> coeff_;
    std::shared_ptr<const Matrix>  nr_coeff_;
    std::shared_ptr<const ZMatrix> hcore_;
    std::shared_ptr<const ZMatrix> overlap_;

    void print_header() const;
    void print_iteration(const int iter, const std::vector<double>& energy, const double error, const double time) const;

    void init();
    virtual void init_mat1e() = 0;
    virtual void init_coeff() = 0;

    // hides some outputs
    mutable std::shared_ptr<Muffle> muffle_;

    std::shared_ptr<ZHarrison> fci_;
    // Fock matrix with active 1RDM
    std::shared_ptr<ZMatrix> compute_active_fock(const ZMatView acoeff, std::shared_ptr<const ZMatrix> rdm1, const bool coulomb_only = false) const;

    // Canonicalize the orbitals
    std::tuple<std::shared_ptr<const ZCoeff_Block>,VectorB,VectorB,VectorB,VectorB> semi_canonical_orb(const bool kramers) const;
    VectorB eig_;
    VectorB eigB_;
    VectorB occup_;
    VectorB occupB_;

    // energy
    std::vector<double> energy_;

    // internal functions
    virtual void impose_symmetry(std::shared_ptr<ZMatrix> o) const = 0;
    virtual void impose_symmetry(std::shared_ptr<ZRotFile> o) const = 0;

    // used to zero out elements for positronic-electronic rotations
    void zero_positronic_elements(std::shared_ptr<ZRotFile> rot);

    std::shared_ptr<ZCoeff_Kramers> nonrel_to_relcoeff(std::shared_ptr<const Matrix> nr_coeff) const;
    std::shared_ptr<const Reference> conv_to_ref_(const bool kramers) const;

  public:
    ZCASSCF(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr);

    virtual void compute() override = 0;

    virtual std::shared_ptr<const Reference> conv_to_ref() const override = 0;

    // functions to retrieve protected members
    int nocc() const { return nocc_; }
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nvirt() const { return nvirt_; }
    int nvirtnr() const { return nvirtnr_; }
    int nbasis() const { return nbasis_; }
    int nstate() const { return nstate_; }
    int charge() const { return charge_; }
    int max_iter() const { return max_iter_; }
    int max_micro_iter() const { return max_micro_iter_; }
    double thresh() const { return thresh_; }
    double thresh_micro() const { return thresh_micro_; }
    double thresh_overlap() const { return thresh_overlap_; }
    double energy(const int i) const { return energy_.at(i); }
    std::vector<double> energy() const { return energy_; }

    std::shared_ptr<const ZHarrison> fci() const { return fci_; }
    std::shared_ptr<const ZCoeff_Block> coeff() const { return coeff_; }
};

}

#endif
