//
// Newint - Parallel electron correlation program.
// Filename: casscf.h
// Copyright (C) 2011 Toru Shiozaki
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
#include <src/wfn/rdm.h>
#include <src/casscf/rotfile.h>

class CASSCF {

  protected:
    // input
    std::multimap<std::string, std::string> idata_; 
    const std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<Reference> ref_;

    // some internal information
    int nocc_; // sum of nact_ + nclosed_
    int nclosed_;
    int nact_;
    int nvirt_;
    // number of MO orbitals. TODO rename to norb. "nbasis" is confusing. 
    int nbasis_;
    int nstate_;
    int max_iter_;
    int max_micro_iter_;
    double thresh_;
    double thresh_micro_;

    std::vector<double> occup_;
    std::shared_ptr<Coeff> coeff_natorb_;

    std::shared_ptr<FCI> fci_;
    void print_header() const;
    void print_iteration(int iter, int miter, int tcount, const std::vector<double> energy, const double error, const double time) const;
    void common_init();

    void mute_stdcout();
    void resume_stdcout();

    std::shared_ptr<Matrix1e> ao_rdm1(std::shared_ptr<RDM<1> > rdm1, const bool inactive_only = false) const;

    const std::shared_ptr<Fock<1> > hcore_;
    void one_body_operators(std::shared_ptr<Matrix1e>&, std::shared_ptr<QFile>&, std::shared_ptr<QFile>&, std::shared_ptr<QFile>&,
                            std::shared_ptr<RotFile>&, const bool superci=true) const;

    std::shared_ptr<const Coeff> update_coeff(const std::shared_ptr<const Coeff>, std::vector<double>) const;
    std::vector<double> form_natural_orbs();

    // energy
    double energy_;

  public:
    CASSCF(const std::multimap<std::string, std::string> idat, const std::shared_ptr<const Geometry> geom);
    virtual ~CASSCF();

    virtual void compute() { assert(false); };

    std::shared_ptr<Reference> ref() { return ref_; };
    std::shared_ptr<const Reference> ref() const { return ref_; };
    virtual std::shared_ptr<const Reference> conv_to_ref() const;

    double energy() const { return energy_; }; 

};

#endif
