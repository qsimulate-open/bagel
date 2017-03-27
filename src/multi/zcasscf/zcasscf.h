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

    // enforce time-reversal symmetry
    bool tsymm_;

    double thresh_;
    double thresh_micro_;

    int nstate_;

    int max_iter_;
    int max_micro_iter_;

    std::shared_ptr<const RelCoeff_Block> coeff_;
    std::shared_ptr<const Matrix>  nr_coeff_;
    std::shared_ptr<const ZMatrix> hcore_;
    std::shared_ptr<const ZMatrix> overlap_;

    void print_header() const;
    void print_iteration(const int iter, const std::vector<double>& energy, const double error, const double time) const;

    void init();

    // hides some outputs
    mutable std::shared_ptr<Muffle> muffle_;

    std::shared_ptr<ZHarrison> fci_;
    // Fock matrix with active 1RDM
    std::shared_ptr<ZMatrix> compute_active_fock(const ZMatView acoeff, std::shared_ptr<const ZMatrix> rdm1, const bool coulomb_only = false) const;

    // energy
    std::vector<double> energy_;

    // internal functions
    // force time-reversal symmetry for a zmatrix with given number of virtual orbitals
    void kramers_adapt(std::shared_ptr<ZMatrix> o, const int nvirt) const;
    // used to zero out elements for positronic-electronic rotations
    void zero_positronic_elements(std::shared_ptr<ZRotFile> rot);

    std::shared_ptr<RelCoeff_Kramers> nonrel_to_relcoeff(std::shared_ptr<const Matrix> nr_coeff) const;

  public:
    ZCASSCF(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref = nullptr);

    virtual void compute() override = 0;

    std::shared_ptr<const Reference> conv_to_ref() const override;

    std::shared_ptr<const RelCoeff_Block> update_coeff(std::shared_ptr<const RelCoeff_Block> cold, std::shared_ptr<const ZMatrix> natorb) const;
    // kramers adapt for RotFile is a static function!
    static void kramers_adapt(std::shared_ptr<ZRotFile> o, const int nclosed, const int nact, const int nvirt);
    // print natural orbital occupation numbers
    void print_natocc(const VectorB& ocup) const;

    // functions to retrieve protected members
    int nocc() const { return nocc_; }
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nvirt() const { return nvirt_; }
    int nvirtnr() const { return nvirtnr_; }
    int nbasis() const { return nbasis_; }
    int nstate() const { return nstate_; }
    int max_iter() const { return max_iter_; }
    int max_micro_iter() const { return max_micro_iter_; }
    double thresh() const { return thresh_; }
    double thresh_micro() const { return thresh_micro_; }
    bool tsymm() const { return tsymm_; }

    std::shared_ptr<const ZHarrison> fci() const { return fci_; }
};

}

#endif
