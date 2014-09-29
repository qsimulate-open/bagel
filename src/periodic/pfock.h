//
// BAGEL - Parallel electron correlation program.
// Filename: pfock.h
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


#ifndef __SRC_PERIODIC_PFOCK_H
#define __SRC_PERIODIC_PFOCK_H

#include <src/periodic/lattice.h>
#include <src/periodic/kdata.h>
#include <src/periodic/pdata.h>

namespace bagel {

class PFock {
  protected:

    std::shared_ptr<PData> pdata_;
    std::shared_ptr<KData> kdata_;

    std::shared_ptr<const Lattice> lattice_;
    std::shared_ptr<const PFock> previous_;
    std::shared_ptr<const ZMatrix> pcoeff_;

    int nblock_, blocksize_;

    // Fourier transform
    void ft();
    // Inverse fourier transform
    void ift();

  public:
    PFock(std::shared_ptr<const Lattice> lattice, std::shared_ptr<const PFock> previous, std::shared_ptr<const ZMatrix> pcoeff);
    ~PFock() { }

    const std::shared_ptr<const PData> pdata() const { return pdata_; }
    const std::shared_ptr<const KData> kdata() const { return kdata_; }
};

}

#endif

