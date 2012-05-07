//
// Newint - Parallel electron correlation program.
// Filename: superci.h
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


#ifndef __NEWINT_CASSCF_SUPERCI_H
#define __NEWINT_CASSCF_SUPERCI_H

#include <memory>
#include <string>
#include <map>
#include <src/scf/scf.h>
#include <src/casscf/casscf.h>
#include <src/wfn/rdm.h>

class SuperCI : public CASSCF {

  protected:
    // DIIS will be used after some macro iteration 
    int diis_start_;

    void common_init() {
      std::cout << "    * Using the Super CI algorithm as noted in Roos (1980) IJQC" << std::endl;
      diis_start_ = read_input<int>(idata_, "diis_start", 5);
      std::cout << "    * DIIS will be used after " << diis_start_ << " macro iteration" << std::endl << std::endl;
    };

    void grad_vc(const std::shared_ptr<Matrix1e> fock, std::shared_ptr<RotFile> sigma);
    void grad_va(const std::shared_ptr<QFile> fact, std::shared_ptr<RotFile> sigma);
    void grad_ca(const std::shared_ptr<Matrix1e> fock, const std::shared_ptr<QFile> fact, std::shared_ptr<RotFile> sigma);

    void compute_qxr(double* int1ext, std::shared_ptr<RDM<2> > rdm2, std::shared_ptr<QFile> qxr);

    void sigma_at_at_(const std::shared_ptr<RotFile> cc, std::shared_ptr<RotFile> sigma,
                      const std::shared_ptr<QFile> gaa, const std::shared_ptr<Matrix1e> f);
    void sigma_ai_ai_(const std::shared_ptr<RotFile> cc, std::shared_ptr<RotFile> sigma, const std::shared_ptr<Matrix1e> f);
    void sigma_at_ai_(const std::shared_ptr<RotFile> cc, std::shared_ptr<RotFile> sigma, const std::shared_ptr<QFile> fact);
    void sigma_ai_ti_(const std::shared_ptr<RotFile> cc, std::shared_ptr<RotFile> sigma, const std::shared_ptr<QFile> fact);
    void sigma_ti_ti_(const std::shared_ptr<RotFile> cc, std::shared_ptr<RotFile> sigma,
                      const std::shared_ptr<QFile> gaa, const std::shared_ptr<Matrix1e> f, const std::shared_ptr<QFile> factp);

    std::shared_ptr<RotFile> const_denom(const std::shared_ptr<QFile> gaa, const std::shared_ptr<QFile> factp,
                                         const std::shared_ptr<Matrix1e> f) const;

    void update_orbitals(std::shared_ptr<RotFile> rot);
    std::shared_ptr<Matrix1e> tailor_rotation(const std::shared_ptr<Matrix1e> seed);

  public:
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref)
      : CASSCF(idat, geom, ref) { common_init(); };
    ~SuperCI() {};

    void compute();

};

#endif
