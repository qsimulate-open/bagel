//
// BAGEL - Parallel electron correlation program.
// Filename: ras/determinants.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ras/determinants.h>

using namespace std;
using namespace bagel;

RASDeterminants::RASDeterminants(const int norb1, const int norb2, const int norb3, const int nelea, const int neleb,
    const int max_holes, const int max_particles, const bool mute) :
  ras_{ norb1, norb2, norb3 }, norb_(norb1 + norb2 + norb3), nelea_(nelea), neleb_(neleb), max_holes_(max_holes), max_particles_(max_particles)
{
  // Construct spaces and with them, a list of strings
  alphaspaces_.reserve( (max_holes_+1) * (max_particles_+1) );
  betaspaces_.reserve( (max_holes_+1) * (max_particles_+1) );

  if (!mute) cout << " o Constructing all possible strings with up to " << max_holes_ << "holes and " << max_particles_ << "particles" << endl;
  for (int nholes = 0; nholes <= max_holes_; ++nholes) {
    const int nele1 = norb1 - nholes;
    for (int nparticles = 0; nparticles <= max_particles_; ++nparticles) {
      const int nele3 = nparticles;
      const int nele2a = nelea_ - (nele1 + nele2);
      const int nele2b = neleb_ - (nele1 + nele2);

      alphaspaces_.push_back( make_shared<const StringSpace(nele1, norb1, nele2a, norb2, nele3, norb3, stringa_.size()) );
      stringa_.insert(stringa_.end(), alphaspaces_.back()->strings().begin(), alphaspaces_.back()->strings().end());

      betaspaces_.push_back( make_shared<const StringSpace(nele1, norb1, nele2b, norb2, nele3, norb3, stringb_.size()) );
      stringb_.insert(stringb_.end(), betaspaces_.back()->strings().begin(), betaspaces_.back()->strings().end());
    }
  }
  if (!mute) cout << "   - alpha strings: " << stringa_.size() << endl;
  if (!mute) cout << "   - beta strings: " << stringb_.size() << endl << endl;

  if (!mute) cout << " o Constructing alpha and beta displacement lists" << endl;
  construct_phis_<0>(stringa_, phia_);
  if (!mute) cout << "   - alpha lists: " << phia_.size() * phia_.front().size() << endl;
  construct_phis_<1>(stringb_, phib_);
  if (!mute) cout << "   - beta lists: " << phib_.size() * phib_.front().size() << endl;
}

#endif
