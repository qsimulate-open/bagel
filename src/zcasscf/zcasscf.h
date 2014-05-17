//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_ZCASSCF_ZCASSCF_H
#define __SRC_ZCASSCF_ZCASSCF_H

#include <src/zfci/zharrison.h>
#include <src/casscf/rotfile.h>
#include <src/wfn/method.h>
#include <src/rel/reloverlap.h>
#include <src/math/bfgs.h>
#include <src/math/step_restrict_bfgs.h>

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

    double thresh_;
    double thresh_micro_;
    std::complex<double> level_shift_;

    int nstate_;

    int max_iter_;
    int max_micro_iter_;

    std::shared_ptr<const ZMatrix> coeff_;
    std::shared_ptr<const Matrix>  nr_coeff_;
    std::shared_ptr<const ZMatrix> hcore_;
    std::shared_ptr<const RelOverlap> overlap_;
    std::vector<double> occup_;

    void print_header() const;
    void print_iteration(int iter, int miter, int tcount, const std::vector<double> energy, const double error, const double time) const;

    void init();
    void init_kramers_coeff();

    void mute_stdcout(const bool fci) const;
    void resume_stdcout() const;

    std::shared_ptr<ZHarrison> fci_;
    std::shared_ptr<const ZMatrix> active_fock(std::shared_ptr<const ZMatrix>, const bool with_hcore = false) const;
    std::shared_ptr<const ZMatrix> transform_rdm1() const;

    // energy
    std::vector<double> energy_;
    double micro_energy_;

    // internal function
    void kramers_adapt(std::shared_ptr<ZRotFile> o, const int nvirt) const;
    void kramers_adapt(std::shared_ptr<ZMatrix> o, const int nvirt) const;

    void zero_positronic_elements(std::shared_ptr<ZRotFile> rot);

    std::shared_ptr<ZMatrix> nonrel_to_relcoeff(const bool stripes = true) const;

  public:
    ZCASSCF(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref = nullptr);

    virtual void compute() override = 0;

    // TODO
    std::shared_ptr<const Reference> conv_to_ref() const override { return nullptr; }

    std::shared_ptr<ZMatrix> make_natural_orbitals(std::shared_ptr<const ZMatrix> rdm1);
    std::shared_ptr<const ZMatrix> natorb_rdm1_transform(const std::shared_ptr<ZMatrix> coeff, std::shared_ptr<const ZMatrix> rdm1) const;
    std::shared_ptr<const ZMatrix> natorb_rdm2_transform(const std::shared_ptr<ZMatrix> coeff, std::shared_ptr<const ZMatrix> rdm1) const;
    std::shared_ptr<const ZMatrix> update_coeff(std::shared_ptr<const ZMatrix> cold, std::shared_ptr<const ZMatrix> natorb) const;
    std::shared_ptr<const ZMatrix> update_qvec(std::shared_ptr<const ZMatrix> qold, std::shared_ptr<const ZMatrix> natorb) const;
    std::shared_ptr<const ZMatrix> semi_canonical_orb();

    // functions to retrieve protected members
    int nocc() const { return nocc_; }
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nvirt() const { return nvirt_; }
    int nbasis() const { return nbasis_; }
    int nstate() const { return nstate_; }
    int max_iter() const { return max_iter_; }
    int max_micro_iter() const { return max_micro_iter_; }
    double thresh() const { return thresh_; }
    double thresh_micro() const { return thresh_micro_; }

};

static const double zoccup_thresh = 1.0e-10;

}

#endif
