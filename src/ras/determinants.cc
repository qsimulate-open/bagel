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

#include <iostream>
#include <iomanip>
#include <cassert>

#include <src/util/combination.hpp>
#include <src/ras/determinants.h>

using namespace std;
using namespace bagel;

RASDeterminants::RASDeterminants(const int norb1, const int norb2, const int norb3, const int nelea, const int neleb,
    const int max_holes, const int max_particles, const bool mute) :
  ras_{{norb1, norb2, norb3}}, norb_(norb1 + norb2 + norb3), nelea_(nelea), neleb_(neleb), max_holes_(max_holes), max_particles_(max_particles),
    lenholes_( ((max_holes_+1)*(max_holes_+2))/2 ), lenparts_( ((max_particles_+1)*(max_particles_+2))/2 )
{
  if ( nelea < 0 || neleb < 0) throw runtime_error("nele < 0");

  // check that large__ is big enough
  if (max_particles_ >= large__) throw logic_error("Inappropriate value for \"large__\". Must be greater than max_particles");

  // Construct spaces and with them, a list of strings
  if (!mute) cout << " o Restricted Active Spaces:" << endl;
  if (!mute) cout << "   - RAS1 -> " << ras_[0] << endl;
  if (!mute) cout << "   - RAS2 -> " << ras_[1] << endl;
  if (!mute) cout << "   - RAS3 -> " << ras_[2] << endl << endl;

  if (!mute) cout << " o Constructing all possible strings with up to " << max_holes_ << " holes and " << max_particles_ << " particles" << endl;
  for (int nholes = 0; nholes <= max_holes_; ++nholes) {
    const int nele1 = norb1 - nholes;
    for (int nparticles = 0; nparticles <= max_particles_; ++nparticles) {
      const int nele3 = nparticles;
      const int nele2a = nelea_ - (nele1 + nele3);
      const int nele2b = neleb_ - (nele1 + nele3);

      if ( (nele1 >= 0) && (nele3 <= norb3) ) {
        if ( (nele2a >= 0) && (nele2a <= norb2) ) {
          auto sp = make_shared<const StringSpace>(nele1, norb1, nele2a, norb2, nele3, norb3, stringa_.size());
          alphaspaces_.emplace(nparticles + nholes * large__, sp);
          stringa_.insert(stringa_.end(), sp->strings().begin(), sp->strings().end());
        }

        if ( (nele2b >= 0) && (nele2b <= norb2) ) {
          auto sp = make_shared<const StringSpace>(nele1, norb1, nele2b, norb2, nele3, norb3, stringb_.size());
          betaspaces_.emplace(nparticles + nholes * large__, sp);
          stringb_.insert(stringb_.end(), sp->strings().begin(), sp->strings().end());
        }
      }
    }
  }
  if (!mute) cout << "   - alpha strings: " << stringa_.size() << endl;
  if (!mute) cout << "   - beta strings: " << stringb_.size() << endl << endl;

  if (!mute) cout << " o Constructing alpha and beta displacement lists" << endl;
  construct_phis_<0>(alphaspaces_, phia_, phia_ij_);
  if (!mute) cout << "   - alpha lists: " << accumulate(phia_.begin(), phia_.end(), 0, [] (const int init, vector<RAS::DMap> plist) { return init + plist.size(); }) << endl;
  construct_phis_<1>(betaspaces_, phib_, phib_ij_);
  if (!mute) cout << "   - beta lists: " << accumulate(phib_.begin(), phib_.end(), 0, [] (const int init, vector<RAS::DMap> plist) { return init + plist.size(); }) << endl;

  if (!mute) cout << " o Constructing pairs of allowed string spaces" << endl;
  size_ = 0;
  stringpairs_.reserve( lenholes_ * lenparts_ );
  for (int nholes = 0; nholes <= max_holes_; ++nholes) {
    for (int nha = nholes; nha >= 0; --nha) {
      const int nhb = nholes - nha;

      for (int npart = 0; npart <= max_particles_; ++npart) {
        for (int npa = npart; npa >= 0; --npa) {
          const int npb = npart - npa;

          shared_ptr<const StringSpace> sa = space<0>(nha, npa);
          shared_ptr<const StringSpace> sb = space<1>(nhb, npb);

          stringpairs_.emplace_back(sa, sb);

          if ( sa && sb ) size_ += sa->size() * sb->size();
        }
      }
    }
  }
  if (!mute) cout << "   - size of restricted space: " << size_ << endl;
}

pair<vector<tuple<bitset<nbit__>, bitset<nbit__>, int>>, double> RASDeterminants::spin_adapt(const int spin, const bitset<nbit__> alpha, const bitset<nbit__> beta) const {
  vector<tuple<bitset<nbit__>, bitset<nbit__>, int>> out;

  // bit pattern for doubly occupied orbitals
  bitset<nbit__> common = (alpha & beta);

  bitset<nbit__> alpha_without_common = alpha ^ common;
  bitset<nbit__> beta_without_common = beta ^ common;

  // alpha pattern without highest spin orbitals
  vector<int> salpha_array = bit_to_numbers(alpha_without_common);
  vector<int> ualpha_array;
  if (salpha_array.size() < spin) throw logic_error("Something is wrong? Determinants::spin_adapt");
  for (int i = 0; i != spin; ++i) {
    ualpha_array.push_back(salpha_array.back());
    salpha_array.pop_back();
  }
  bitset<nbit__> salpha = numbers_to_bit(salpha_array);
  bitset<nbit__> ualpha = numbers_to_bit(ualpha_array);
  bitset<nbit__> common_plus_alpha(common.to_ulong() + ualpha.to_ulong());

  // number of unpaired alpha orbitals (minus Ms)
  const int nalpha = salpha.count();

  // a vector of number that specify open orbitals
  vector<int> open = bit_to_numbers(salpha^beta_without_common);
  assert((salpha^beta_without_common) == (salpha|beta_without_common));

  // take a linear combination to make a vector singlet coupled.
  // TODO for the time being, we just leave Ms highest orbitals and singlet-couple other orbitals
  assert(nalpha*2 == open.size());
  int icnt = 0;
  do {
    bitset<nbit__> ialpha = common_plus_alpha;
    bitset<nbit__> ibeta = common;
    for (int i =0; i!=nalpha; ++i) ialpha.flip(open[i]);
    for (int i=nalpha; i!=open.size(); ++i) ibeta.flip(open[i]);

    // sign is always compensated by moving alpha to the left and beta to the right
    // our convention is (aaaa)(bbbb)|0> due to the alpha-beta string algorithm
    const double sign = 1.0;

    out.push_back(make_tuple(ibeta, ialpha, sign));
    ++icnt;
  } while (boost::next_combination(open.begin(), open.begin()+nalpha, open.end()));

  // scale to make the vector normalized
  const double factor = 1.0/sqrt(static_cast<double>(icnt));
  return make_pair(out, factor);
}
