//
// BAGEL - Parallel electron correlation program.
// Filename: space_base.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
// Modified by: Michael Caldwell <caldwell@u.northwestern.edu>
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


#ifndef __SRC_FCI_SPACE_BASE_H
#define __SRC_FCI_SPACE_BASE_H

#include <src/fci/determinants.h>

namespace bagel {

// implements a space that contains multiple determinants made by modification of a reference
class Space_base {
  protected:
    // assuming that the number of active orbitals are the same in alpha and beta.
    const int norb_;

    const int nelea_; // reference number of alpha electrons
    const int neleb_; // reference number of beta electrons

    const bool mute_;

    int key_(const int a, const int b) const { return ( a*large__ + b ); }

    std::map<int, std::shared_ptr<Determinants>> detmap_; // For now, all access should be through Determinants objects

    virtual void common_init() = 0;

  public:
    Space_base(std::shared_ptr<const Determinants> det_, const bool mute = true) :
      norb_(det_->norb()), nelea_(det_->nelea()), neleb_(det_->neleb()), mute_(mute) { }
    Space_base(const int norb, const int nelea, const int neleb, const bool mute = true) :
      norb_(norb), nelea_(nelea), neleb_(neleb), mute_(mute) { }

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    std::shared_ptr<Determinants> basedet() { return finddet(0, 0); };
    std::shared_ptr<const Determinants> basedet() const { return finddet(0, 0); };
    // Caution: This function does not check to make sure i,j is valid
    std::shared_ptr<Determinants> finddet(const int i, const int j) { auto idet = detmap_.find(key_(i,j)); return idet->second; };
    std::shared_ptr<const Determinants> finddet(const int i, const int j) const { auto idet = detmap_.find(key_(i,j)); return idet->second; };

    const std::map<int, std::shared_ptr<Determinants>>& detmap() const { return detmap_; }

    int key(std::shared_ptr<const Determinants> det) const { return key_(det->nelea() - nelea_, det->neleb() - neleb_); }
};

}

#endif
