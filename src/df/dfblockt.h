//
// BAGEL - Parallel electron correlation program.
// Filename: dfblockt.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __SRC_DF_DFBLOCKT_H
#define __SRC_DF_DFBLOCKT_H

#include <map>
#include <src/df/dfblock.h>

namespace bagel {

/*
    DFBlockT is a slice of 3-index DF integrals. Distributed by the second and third index.
    Note that the date is tranposed!
*/

class DFBlockT {
  protected:
    std::unique_ptr<double[]> data_;

    // first dimension is naux_ (global)
    const size_t naux_;
    // second and third dimension
    const size_t bstart_;
    const size_t bsize_;

    // map between an index number and b1 and b2 number
    const std::map<size_t, std::pair<size_t, size_t> > index_;

  public:
    DFBlockT(std::shared_ptr<const DFBlock> in, const size_t naux, const size_t bstart, const size_t bsize, const std::map<size_t, std::pair<size_t, size_t> >& index);

}; 
    
}

#endif
