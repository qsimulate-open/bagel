//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcassecond.h
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

#ifndef __SRC_ZCASSCF_ZCASSECOND_H
#define __SRC_ZCASSCF_ZCASSECOND_H

#include <src/multi/zcasscf/zcasscf.h>
#include <src/df/reldfhalf.h>

namespace bagel {

class ZCASSecond_base : public ZCASSCF {
  protected:
    // convergence threshold for micro iteration relative to stepsize
    double thresh_microstep_;

    // print out converged orbitals in the DIRAC's DFPCMO format
    bool dfpcmo_;

    // compute orbital gradient
    std::shared_ptr<ZRotFile> compute_gradient(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr) const;

    // diagonal Hessian
    std::shared_ptr<ZRotFile> compute_denom(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock,
                                            std::shared_ptr<const ZMatrix> qxr, std::shared_ptr<const ZMatrix> rdm1) const;
    // compute H*t (Hessian times trial vector)
    std::shared_ptr<ZRotFile>
      compute_hess_trial(std::shared_ptr<const ZRotFile> trot, std::list<std::shared_ptr<const RelDFHalf>> halfc, std::list<std::shared_ptr<const RelDFHalf>> halfac,
                         std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr) const;
    // apply denominator in microiterations
    std::shared_ptr<ZRotFile> apply_denom(std::shared_ptr<const ZRotFile> grad, std::shared_ptr<const ZRotFile> denom, const double shift, const double scale) const;

    // init functions
    virtual void init_mat1e() override = 0;
    virtual void init_coeff() override = 0;

    virtual void impose_symmetry(std::shared_ptr<ZMatrix>) const override = 0;
    virtual void impose_symmetry(std::shared_ptr<ZRotFile>) const override = 0;

    virtual void trans_natorb() = 0;
    virtual bool kramers() const = 0;

  protected:
    ZCASSecond_base(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref);

  public:
    void compute() override;
    virtual std::shared_ptr<const Reference> conv_to_ref() const override = 0;
};


class ZCASSecond : public ZCASSecond_base {
  protected:
    virtual void init_mat1e() override final;
    virtual void init_coeff() override final;
    virtual void impose_symmetry(std::shared_ptr<ZMatrix>) const override final;
    virtual void impose_symmetry(std::shared_ptr<ZRotFile>) const override final;
    virtual void trans_natorb() override final;
    virtual bool kramers() const override { return true; }

  public:
    ZCASSecond(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr);
    std::shared_ptr<const Reference> conv_to_ref() const override { return conv_to_ref_(true); }
};


class ZCASSecond_London : public ZCASSecond_base {
  protected:
    virtual void init_mat1e() override final;
    virtual void init_coeff() override final;
    virtual void impose_symmetry(std::shared_ptr<ZMatrix>) const override final { /*do nothing*/ }
    virtual void impose_symmetry(std::shared_ptr<ZRotFile>) const override final { /*do nothing*/ }
    virtual void trans_natorb() override final;
    virtual bool kramers() const override { return false; }

  public:
    ZCASSecond_London(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr);
    std::shared_ptr<const Reference> conv_to_ref() const override { return conv_to_ref_(false); }
};

}

#endif
