//
// Newint - Parallel electron correlation program.
// Filename: scf_base.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#ifndef __scf_scf_base_h
#define __scf_scf_base_h

#include <memory>
#include <string>
#include <map>
#include <src/scf/geometry.h>
#include <src/scf/overlap.h>
#include <src/scf/hcore.h>
#include <src/scf/tildex.h>
#include <src/scf/fock.h>
#include <src/scf/coeff.h>
#include <src/wfn/reference.h>

class SCF_base {
  protected:
    std::multimap<std::string, std::string> idata_;
    const std::shared_ptr<Geometry> geom_;
    const std::shared_ptr<Overlap> overlap_;
    const std::shared_ptr<Hcore> hcore_;
    std::shared_ptr<TildeX> tildex_;
    std::shared_ptr<Matrix1e> aodensity_;
    std::shared_ptr<Coeff> coeff_;

    int max_iter_;
    int diis_start_;
    double thresh_overlap_;
    double thresh_scf_;
    bool density_change_;

    std::vector<double> schwarz_;
    void init_schwarz();

    std::unique_ptr<double[]> eig_;

    int nocc_;

  public:
    SCF_base(std::multimap<std::string, std::string>& idata_, const std::shared_ptr<Geometry>);
    ~SCF_base() {};

    virtual void compute() = 0;

    const std::shared_ptr<Geometry> geom() { return geom_; };
    const std::shared_ptr<Matrix1e> aodensity() { return aodensity_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    void set_coeff(const std::shared_ptr<Coeff> o) { coeff_ = o; };
    const std::shared_ptr<Hcore> hcore() { return hcore_; };
    const std::vector<double>& schwarz() const { return schwarz_; };

    int nocc() const { return nocc_; };

    std::shared_ptr<Reference> conv_to_ref() {
      std::shared_ptr<Reference> out(new Reference(geom_, coeff(), hcore(), schwarz(), nocc(), 0, geom_->nbasis()-nocc()));
      std::vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
      out->set_eig(e);
      return out;
    };

    double* eig() { return eig_.get(); };
};

#endif
