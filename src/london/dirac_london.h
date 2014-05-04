//
// BAGEL - Parallel electron correlation program.
// Filename: dirac_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_LONDON_DIRAC_LONDON_H
#define __BAGEL_SRC_LONDON_DIRAC_LONDON_H

#include <string>
#include <map>
#include <src/wfn/method.h>
#include <src/wfn/geometry_london.h>
#include <src/london/relhcore_london.h>
#include <src/london/reloverlap_london.h>
#include <src/london/reldfhalf_london.h>

namespace bagel {

class Dirac_London : public Method {
  protected:

    int max_iter_;
    int diis_start_;
    double thresh_scf_;
    double thresh_overlap_;
    double energy_;
    std::unique_ptr<double[]> eig_;
    int ncharge_;
    int nele_;
    int nneg_;

    bool gaunt_;
    bool breit_;

    // for Fock build
    bool robust_;

    std::shared_ptr<const RelHcore_London> hcore_;
    std::shared_ptr<const RelOverlap_London> overlap_;
    std::shared_ptr<const ZMatrix> s12_;

    std::shared_ptr<const ZMatrix> coeff_;

    void common_init(const std::shared_ptr<const PTree>);
    std::shared_ptr<const DistZMatrix> initial_guess(const std::shared_ptr<const DistZMatrix> s12, const std::shared_ptr<const DistZMatrix> hcore) const;

    // if gradient is requested, half-transformed integrals will be reused
    bool do_grad_;
    std::list<std::shared_ptr<RelDFHalf_London>> half_;

  public:
    Dirac_London() { }
    Dirac_London(const std::shared_ptr<const PTree> idata_, const std::shared_ptr<const Geometry_London> geom, const std::shared_ptr<const Reference> re = nullptr);

    ~Dirac_London() { geom_->discard_relativistic(); }

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;

    void print_eig() const;
    double energy() const { return energy_; }

    std::shared_ptr<const Geometry_London> geom() const { return cgeom_; }

    std::list<std::shared_ptr<RelDFHalf_London>> half() const { return half_; }
    void discard_half() { half_.clear(); }

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Dirac_London)

#endif
