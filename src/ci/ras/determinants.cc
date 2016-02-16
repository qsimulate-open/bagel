//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/determinants.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <iostream>
#include <iomanip>
#include <cassert>

#include <src/util/combination.hpp>
#include <src/ci/ras/determinants.h>

using namespace std;
using namespace bagel;

RASDeterminants::RASDeterminants(const int norb1, const int norb2, const int norb3, const int nelea, const int neleb,
                                 const int max_holes, const int max_particles, const bool mute)
 : ras_{{norb1, norb2, norb3}}, max_holes_(max_holes), max_particles_(max_particles) {

  if ( nelea < 0 || neleb < 0) throw runtime_error("nele < 0");

  // check that nbit__ is big enough
  if ( max_particles_ >= nbit__ || std::max(nelea,neleb) >= nbit__)
    throw logic_error("Inappropriate value for \"nbit__\". Must be greater than nele AND max_particles");

  // Construct spaces and with them, a list of strings
  if (!mute) cout << " o Restricted Active Spaces:" << endl;
  if (!mute) cout << "   - RAS1 -> " << ras_[0] << endl;
  if (!mute) cout << "   - RAS2 -> " << ras_[1] << endl;
  if (!mute) cout << "   - RAS3 -> " << ras_[2] << endl << endl;

  if (!mute) cout << " o Constructing all possible strings with up to " << max_holes_ << " holes and " << max_particles_ << " particles" << endl;
  {
    list<shared_ptr<const RASString>> alpha;
    list<shared_ptr<const RASString>> beta;

    for (int nholes = 0; nholes <= max_holes_; ++nholes) {
      const int nele1 = norb1 - nholes;
      for (int nparticles = 0; nparticles <= max_particles_; ++nparticles) {
        const int nele3 = nparticles;
        const int nele2a = nelea - (nele1 + nele3);
        const int nele2b = neleb - (nele1 + nele3);

        if (nele1 >= 0 && nele3 <= norb3) {
          if (nele2a >= 0 && nele2a <= norb2) {
            auto astring = make_shared<RASString>(nele1, norb1, nele2a, norb2, nele3, norb3);
            if (astring->size() > 0)
              alpha.push_back(astring);
          }

          if (nele2b >= 0 && nele2b <= norb2) {
            auto bstring = make_shared<RASString>(nele1, norb1, nele2b, norb2, nele3, norb3);
            if (bstring->size() > 0)
              beta.push_back(bstring);
          }
        }
      }
    }
    if (!alpha.empty())
      alphaspaces_ = make_shared<CIStringSet<RASString>>(alpha);
    if (!beta.empty())
      betaspaces_ = make_shared<CIStringSet<RASString>>(beta);
  }

  size_ = 0;

  if (alphaspaces_ && betaspaces_) {
    if (alphaspaces_->size()*betaspaces_->size() > 0) {
      if (!mute) cout << "   - alpha strings: " << alphaspaces_->size() << endl;
      if (!mute) cout << "   - beta strings: "  <<  betaspaces_->size() << endl << endl;

      if (!mute) cout << " o Constructing alpha and beta displacement lists" << endl;
      construct_phis_<0>(alphaspaces_, phia_, phia_ij_);
      if (!mute) cout << "   - alpha lists: " << phia_->size() << endl;
      construct_phis_<1>(betaspaces_, phib_, phib_ij_);
      if (!mute) cout << "   - beta lists: "  << phib_->size() << endl;

      if (!mute) cout << " o Constructing pairs of allowed string spaces" << endl;

      for (int nholes = 0; nholes <= max_holes_; ++nholes) {
        for (int nha = nholes; nha >= 0; --nha) {
          const int nhb = nholes - nha;
          for (int npart = 0; npart <= max_particles_; ++npart) {
            for (int npa = npart; npa >= 0; --npa) {
              const int npb = npart - npa;

              auto block = make_shared<const CIBlockInfo<RASString>>(space<0>(nha, npa), space<1>(nhb, npb), size_);
              blockinfo_.push_back(block);
              if (!block->empty()) size_ += block->size();
            }
          }
        }
      }
    }
    if (!mute) cout << "   - size of restricted space: " << size_ << endl;
  }
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
  return {out, factor};
}
