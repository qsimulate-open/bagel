//
// Newint - Parallel electron correlation program.
// Filename: casscf.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

// This is a base class for various CASSCF solvers.
// The assumpation is made that the CI and orbital optimizations are done independently.
// This should be a good strategy in large systems.
//

#ifndef __NEWINT_CASSCF_CASSCF_H
#define __NEWINT_CASSCF_CASSCF_H

#include <cassert>
#include <vector>
#include <memory>
#include <src/scf/geometry.h>
#include <src/wfn/reference.h>
#include <src/fci/fci.h>
#include <src/fci/rdm.h>
#include <src/casscf/rotfile.h>

class CASSCF {

  protected:
    // input
    std::multimap<std::string, std::string> idata_; 

    // some internal information
    int nocc_; // sum of nact_ + nclosed_
    int nclosed_;
    int nact_;
    int nvirt_;
    int nbasis_;
    int nstate_;
    int max_iter_;
    int max_micro_iter_;
    double thresh_;
    double thresh_micro_;

    std::vector<double> occup_;
    std::shared_ptr<Coeff> coeff_natorb_;

    std::shared_ptr<Reference> ref_;
    std::shared_ptr<FCI> fci_;
    const std::shared_ptr<Geometry> geom_;
    void print_header() const;
    void common_init();

    void mute_stdcout();
    void resume_stdcout();

    std::shared_ptr<Matrix1e> ao_rdm1(std::shared_ptr<RDM<1> > rdm1, const bool inactive_only = false) const;

    std::shared_ptr<Fock<1> > hcore_;
    void one_body_operators(std::shared_ptr<Matrix1e>&, std::shared_ptr<QFile>&, std::shared_ptr<QFile>&, std::shared_ptr<QFile>&,
                            std::shared_ptr<RotFile>&, const bool superci=true);

    std::shared_ptr<Coeff> update_coeff(const std::shared_ptr<Coeff>, std::vector<double>) const;
    std::vector<double> form_natural_orbs();

  public:
    CASSCF(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref);
    virtual ~CASSCF();

    virtual void compute() { assert(false); };

    std::shared_ptr<Reference> ref() { return ref_; };
    std::shared_ptr<Reference> conv_to_ref() const {
      std::shared_ptr<Reference> out(new Reference(geom_, ref_->coeff(), ref_->hcore(), ref_->schwarz(), nclosed_, nact_, nvirt_));
      return out;
    };

};

#endif
