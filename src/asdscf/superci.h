//
// BAGEL - Parallel electron correlation program.
// Filename: superci.h
// Copyright (C) 2011 Toru Shiozaki
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


#ifndef __BAGEL_ASDSCF_SUPERCI_H
#define __BAGEL_ASDSCF_SUPERCI_H

#include <map>
#include <src/scf/scf.h>
#include <src/asdscf/asdscf.h>
#include <src/wfn/rdm.h>

namespace bagel {


class ASDSuperCI : public ASDSCF {

  protected:
    // DIIS will be used after some macro iteration
    int diis_start_;

    void common_init() {
      std::cout << "    * Using the Super CI algorithm as noted in Roos (1980) IJQC" << std::endl;
//TODO
      diis_start_ = 5; //idata_->get<int>("diis_start", 5);
      std::cout << "    * DIIS will be used after " << diis_start_ << " macro iteration" << std::endl << std::endl;
    }

    void grad_vc(const std::shared_ptr<Matrix> fock, std::shared_ptr<ASDRotFile> sigma);
    void grad_va(const std::shared_ptr<Matrix> fact, std::shared_ptr<ASDRotFile> sigma);
    void grad_ca(const std::shared_ptr<Matrix> fock, const std::shared_ptr<Matrix> fact, std::shared_ptr<ASDRotFile> sigma);


    void update_orbitals(std::shared_ptr<ASDRotFile> rot);
    std::shared_ptr<Matrix> tailor_rotation(const std::shared_ptr<Matrix> seed);

  public:
  //ASDSuperCI(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
    ASDSuperCI(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, 
               std::tuple<std::shared_ptr<RDM<1>>,
                          std::shared_ptr<RDM<2>>> rdms )
      : ASDSCF(idat, geom, ref, rdms) { common_init(); }

    void compute() override;

};

}

#endif
