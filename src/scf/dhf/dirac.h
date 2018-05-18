//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dirac.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_REL_DIRAC_H
#define __SRC_REL_DIRAC_H

#include <src/wfn/method.h>
#include <src/df/reldfhalf.h>

namespace bagel {

class Dirac : public Method {
  protected:

    int max_iter_;
    int diis_start_;
    double thresh_scf_;
    double thresh_overlap_;
    double energy_;
    VectorB eig_;
    int ncharge_;
    int nele_;
    int nneg_;

    bool gaunt_;
    bool breit_;

    bool dfpcmo_;

    // for Fock build
    bool robust_;

    int multipole_print_;
    bool conv_ignore_;

    std::shared_ptr<const ZMatrix> hcore_;
    std::shared_ptr<const ZMatrix> overlap_;
    std::shared_ptr<const ZMatrix> s12_;

    std::shared_ptr<const ZMatrix> coeff_;

    void common_init(const std::shared_ptr<const PTree>);
    std::shared_ptr<const DistZMatrix> initial_guess(const std::shared_ptr<const DistZMatrix> s12, const std::shared_ptr<const DistZMatrix> hcore) const;

    // if gradient is requested, half-transformed integrals will be reused
    bool do_grad_;
    std::list<std::shared_ptr<RelDFHalf>> half_;

  public:
    Dirac(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re = nullptr);

    ~Dirac() { geom_->discard_relativistic(); }

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;

    void print_eig() const;
    double energy() const { return energy_; }

    std::shared_ptr<const Geometry> geom() const { return geom_; }

    std::list<std::shared_ptr<RelDFHalf>> half() const { return half_; }
    void discard_half() { half_.clear(); }

};

}

#endif
