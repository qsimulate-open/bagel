//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casbfgs.h
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


#ifndef __BAGEL_CASSCF_CASBFGS_H
#define __BAGEL_CASSCF_CASBFGS_H

#include <src/multi/casscf/casscf.h>

namespace bagel {

// BFGS driver
class CASBFGS : public Method {
  protected:
    std::vector<double> energy_;

    std::shared_ptr<const Reference> refout_;
    std::shared_ptr<FCI> fci_;
    double rms_grad_;

  public:
    CASBFGS(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
     : Method(idat, geom, ref) { }

    void compute() override;
    std::shared_ptr<const Reference> conv_to_ref() const override { assert(refout_); return refout_; }

    double energy(const int i) const { return energy_[i]; }
    const std::vector<double>& energy() const { return energy_; }

    std::shared_ptr<FCI> fci() { return fci_; }
    double rms_grad() const { return rms_grad_; }
};

class CASBFGS_base : public CASSCF {

  protected:
    void common_init() {
      std::cout << "    * Using the Quasi 2nd-order algorithm as noted in Chaban et al. TCA (1997)" << std::endl << std::endl;
    }

    // compute orbital gradients
    void grad_vc(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<RotFile> sigma) const;
    void grad_va(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> qxr,   std::shared_ptr<RotFile> sigma) const;
    void grad_ca(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<RotFile> sigma) const;

    // compute diagonal denominators
    std::shared_ptr<const RotFile> compute_denom(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr) const;

    CASBFGS_base(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
      : CASSCF(idat, geom, ref) { common_init(); }

    virtual void compute() override = 0;
};


// uses BAGEL's native BFGS
class CASBFGS1 : public CASBFGS_base {
  public:
    CASBFGS1(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
     : CASBFGS_base(idat, geom, ref) { }
    void compute() override;
};


// uses alglib's BFGS
class CASBFGS2 : public CASBFGS_base {
  protected:
    bool only_energy_converged_;
  public:
    CASBFGS2(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
     : CASBFGS_base(idat, geom, ref), only_energy_converged_(false) { }

    void compute() override;
    bool only_energy_converged() const { return only_energy_converged_; }
};

}

#endif
