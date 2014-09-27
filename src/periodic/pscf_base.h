//
// BAGEL - Parallel electron correlation program.
// Filename: pscf_base.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#ifndef __BAGEL_SRC_PERIODIC_PSCF_BASE_H
#define __BAGEL_SRC_PERIODIC_PSCF_BASE_H

#include <src/wfn/method.h>
#include <src/periodic/lattice.h>
#include <src/periodic/phcore.h>
#include <src/periodic/poverlap.h>
#include <src/periodic/pcoeff.h>

namespace bagel {

class PSCF_base : public Method {
  protected:
    std::shared_ptr<const Lattice> lattice_;
    std::shared_ptr<const Data> tildex_;
    std::shared_ptr<const POverlap> overlap_;
    std::shared_ptr<const PHcore> hcore_;
    std::shared_ptr<const PCoeff> coeff_;

    int max_iter_;

    int diis_start_;
    int diis_size_;

    double thresh_overlap_;
    double thresh_scf_;

    VectorB eig_;
    double energy_;

    int nocc_;
    int noccB_;

    const std::string indent = "  ";

  public:
    PSCF_base() { }
    PSCF_base(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom,
              const std::shared_ptr<const Reference> re = nullptr);

    virtual ~PSCF_base() { }

    virtual void compute() override = 0;

    const std::shared_ptr<const PCoeff> coeff() const { return coeff_; }
    void set_coeff(const std::shared_ptr<PCoeff> o) { coeff_ = o; };

    const std::shared_ptr<const PHcore> hcore() const { return hcore_; }

    int nocc() const { return nocc_; }
    int noccB() const { return noccB_; }
    double energy() const { return energy_; }

    double thresh_overlap() const { return thresh_overlap_; }
    double thresh_scf() const { return thresh_scf_; }

    virtual std::shared_ptr<const Reference> conv_to_ref() const override = 0;

    VectorB& eig() { return eig_; }
};

}

#endif

