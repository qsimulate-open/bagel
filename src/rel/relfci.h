//
// BAGEL - Parallel electron correlation program.
// Filename: relfci.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Matthew Kelley matthewkelley2013@u.northwestern.edu
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

#ifndef __BAGEL_RELFCI_FCI_H
#define __BAGEL_RELFCI_FCI_H

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

namespace bagel {

class RelFCI {

  protected:
    // input
    std::multimap<std::string, std::string> idata_; 
    // reference
    std::shared_ptr<const RelReference> ref_;
    // geometry file
    const std::shared_ptr<const Geometry> geom_;
    // max #iteration
    int max_iter_;
    // threshold for variants
    double thresh_;
    double print_thresh_;

  public:
    // this constructor is ugly... to be fixed some day...
    RelFCI(const std::multimap<std::string, std::string>, std::shared_ptr<const RelReference>);
        

    void compute();

};

}

#endif

