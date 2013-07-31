//
// BAGEL - Parallel electron correlation program.
// Filename: space.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_ZFCI_RELSPACE_H
#define __SRC_ZFCI_RELSPACE_H

#include <memory>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <bitset>
#include <cassert>

#include <src/util/constants.h>
#include <src/fci/determinants.h>

namespace bagel {

/************************************************************************************
*     Note: I've been using as many of the member functions of bitset as possible,  *
*        not necessarily because I think it will be faster (I have no idea) but     *
*        because I just want to. I can always change back to faster routines        *
************************************************************************************/

// implements spaces that contain all determinants |PQ> for a given Kramers index -N/2 to N/2
class RelSpace {
  protected:
    // assuming that the number of active orbitals are the same in "alpha" and "beta."
    const int norb_;

    const int nunbar_; // reference number of unbarred electrons
    const int nbar_; // reference number of barred electrons

    const bool mute_;

    std::map<int, std::shared_ptr<Determinants>> detmap_; // For now, all access should be through Determinants objects

    int nelec() { return nunbar_ + nbar_; }
    int kramers(int m) { return (nunbar_ - nbar_ - 2*m)/2; }

  public:
#if 0
    Space(std::shared_ptr<const Determinants>, const int M, const bool compress = false, const bool mute = true);
#endif
    RelSpace(const int norb, const int nunbar, const int nbar, const bool mute = true);
    ~RelSpace() {};

    // static constants
    static const int Unbar = 0;
    static const int Bar = 1;

    std::shared_ptr<Determinants> basedet() { return finddet(0); };

    std::shared_ptr<Determinants> finddet(int m) {
      // TODO only for temporary debugging purposes
     assert(abs(m)<=abs(nelec()/2));
      auto idet = detmap_.find(m); return idet->second; };

  private:
    void common_init();
};

}

#endif
