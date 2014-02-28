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

namespace bagel {

class ZCASSCF : public Method {
  protected:
    int nneg_;
    int nocc_;
    int nclosed_;
    int nact_;
    int nvirt_;
    int nbasis_;

    int charge_;

    bool gaunt_;
    bool breit_;

    double thresh_;
    double thresh_micro_;

    int nstate_;

    int max_iter_;
    int max_micro_iter_;

    std::shared_ptr<const ZMatrix> coeff_;

    void print_header() const;
    void print_iteration(int iter, int miter, int tcount, const std::vector<double> energy, const double error, const double time) const;

    void init();
    void init_kramers_coeff(std::shared_ptr<const ZMatrix> hcore, std::shared_ptr<const RelOverlap> overlap);

    void mute_stdcout() const;
    void resume_stdcout() const;

    std::shared_ptr<ZHarrison> fci_;
    std::shared_ptr<const ZMatrix> active_fock(std::shared_ptr<const ZMatrix>) const;
    std::shared_ptr<const ZMatrix> transform_rdm1() const;

    // energy
    std::vector<double> energy_;

    // internal function
    void grad_vc(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<ZRotFile> sigma) const;
    void grad_va(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> qxr,   std::shared_ptr<const ZMatrix> rdm1, std::shared_ptr<ZRotFile> sigma) const;
    void grad_ca(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr,
                 std::shared_ptr<const ZMatrix> rdm1, std::shared_ptr<ZRotFile> sigma) const;

    std::shared_ptr<const ZRotFile> compute_denom(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock,
                                                  std::shared_ptr<const ZMatrix> qxr, std::shared_ptr<const ZMatrix> rdm1) const;

    void kramers_adapt(std::shared_ptr<ZRotFile> o) const;
    void kramers_adapt(std::shared_ptr<ZMatrix> o) const;

  public:
    ZCASSCF(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom,
            const std::shared_ptr<const Reference> = std::shared_ptr<const Reference>());

    void compute();

    // TODO
    std::shared_ptr<const Reference> conv_to_ref() const override { return std::shared_ptr<const Reference>(); }

  private:
    // TODO debug only. All implemented in zcasscf_debug.cc. Will be removed once everything works.
    void ___debug___orbital_rotation(const bool kramers);
    void ___debug___print_gradient(std::shared_ptr<const ZRotFile> grad, const bool kramers) const;
    void ___debug___compute_hessian(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock) const;
    // returns [x,y] = (xx|yy) (x is an index of coeffa, and y is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_coulomb(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,y] = (xy|yx) (x is an index of coeffa, and y is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_exchange(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,t,u] = (xx|tu) (x is an index of coeffa, and coeffi should be active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_coulomb_active(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,t,u] = (xu|tx) (x is an index of coeffa, and coeffi should be active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_exchange_active(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;

};

}

#endif
