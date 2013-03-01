
// BAGEL - Parallel electron correlation program.
// Filename: dimer_cispace.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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



#ifndef __BAGEL_DIMER_CISPACE_H
#define __BAGEL_DIMER_CISPACE_H

#include <fstream>
#include <string>
#include <memory>
#include <cassert>
#include <tuple>
#include <src/wfn/reference.h>
#include <src/fci/mofile.h>

namespace bagel {

class DimerCISpace {
  protected:
    std::shared_ptr<Dimer> dimer_;

    // These are stored values of the neutral species
    std::pair<int, int> nelea_;
    std::pair<int, int> neleb_;

    std::multimap<std::pair<int,int>, std::shared_ptr<Dvec>> cispaceA_;
    std::multimap<std::pair<int,int>, std::shared_ptr<Dvec>> cispaceB_;

  public:
    // This constructor will build the infrastructure; civecs need to be added later
    DimerCISpace(const std::shared_ptr<Dimer> dimer);

    template<int unit> int nelea() { return (unit == 0 ? nelea_.first : nelea_.second); }
    template<int unit> int neleb() { return (unit == 0 ? neleb_.first : neleb_.second); }

    template<int unit> void insert(std::shared_ptr<const Dvec> civec);
};

template<int unit> void insert(std::shared_ptr<const Dvec> civec) {
  std::shared_ptr<Dvec> new_civec = civec->copy();

  // Reform Determinants object (to make sure it's the format I want)
  std::shared_ptr<Determinants> det(new Determinants(civec->det(), /*compress=*/false, /*mute=*/true));
  new_civec->set_det(det);

  const int qa = det->nelea() - nelea<unit>();
  const int qb = det->neleb() - neleb<unit>();

  auto& space = (unit == 0 ? cispaceA_ : cispaceB_);

  auto icivec = space.find(make_pair(qa+1,qb));
  if (icivec != space.end()) det->link<0>((*icivec)->det());

  auto icivec = space.find(make_pair(qa-1,qb));
  if (icivec != space.end()) det->link<0>((*icivec)->det());

  auto icivec = space.find(make_pair(qa,qb+1));
  if (icivec != space.end()) det->link<1>((*icivec)->det());

  auto icivec = space.find(make_pair(qa,qb-1));
  if (icivec != space.end()) det->link<1>((*icivec)->det());

  space.insert(make_pair(qa,qb), new_civec);
}

}

#endif
