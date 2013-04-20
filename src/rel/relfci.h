//
// BAGEL - Parallel electron correlation program.
// Filename: relfci.h
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2013@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __BAGEL_FCI_RELFCI_H
#define __BAGEL_FCI_RELFCI_H

#include <tuple>
#include <cassert>
#include <iostream>
#include <memory>
#include <bitset>
#include <src/util/constants.h>
#include <src/fci/dvec.h>
#include <src/fci/mofile.h>
#include <src/fci/determinants.h>
#include <src/fci/properties.h>
#include <src/wfn/rdm.h>
#include <src/wfn/reference.h>
#include <src/rel/relreference.h>
#include <src/wfn/ciwfn.h>
#include <src/util/zmatrix.h>

namespace bagel {

class RelFCI {

  protected:
    std::multimap<std::string, std::string> idata_;

    double thresh_;
    double print_thresh_;
    std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const RelReference> relref_;

    int max_iter_;
    int diis_start_;
    double thresh_scf_;
    double energy_;
    int ncharge_;
    int nele_;
    int nneg_;

    bool gaunt_;
    bool breit_;

    std::shared_ptr<const ZMatrix> coeff_;

    void common_init(const std::multimap<std::string, std::string>&);

    std::shared_ptr<const ZMatrix> time_reversal_operator();

  public:
    RelFCI(const std::multimap<std::string, std::string>&, const std::shared_ptr<const Geometry> geom,
             const std::shared_ptr<const RelReference> re);

    void compute();

    void print_eig(const std::unique_ptr<double[]>&);

};

}

#endif

