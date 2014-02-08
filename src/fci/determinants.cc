//
// BAGEL - Parallel electron correlation program.
// Filename: determinants.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <stdexcept>
#include <iomanip>
#include <src/fci/determinants.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Determinants)

using namespace std;
using namespace bagel;

Determinants::Determinants(const int _norb, const int _nelea, const int _neleb, const bool _compress, const bool mute) : compress_(_compress) {

  blockinfo_[0] = make_shared<FCIBlockInfo>(_norb, _nelea, _neleb, mute);

  for (auto& i : blockinfo_) {
    string_bits_a_.insert(string_bits_a_.end(), i.second->string_bits_a().begin(), i.second->string_bits_a().end());
    string_bits_b_.insert(string_bits_b_.end(), i.second->string_bits_b().begin(), i.second->string_bits_b().end());
  }

  if (!mute) cout << "  o single displacement lists (alpha)" << endl;
  const_phis_<0>(string_bits_a(), phia_, phia_uncompressed_);
  if (!mute) cout << "      length: " << setw(13) << accumulate(phia_.begin(), phia_.end(), 0, [](const int init, vector<DetMap>& plist) { return init + plist.size(); }) << endl;
  if (!mute) cout << "  o single displacement lists (beta)" << endl;
  const_phis_<1>(string_bits_b(), phib_, phib_uncompressed_);
  if (!mute) cout << "      length: " << setw(13) << accumulate(phib_.begin(), phib_.end(), 0, [](const int init, vector<DetMap>& plist) { return init + plist.size(); }) << endl;

}


pair<vector<tuple<int, int, int>>, double> Determinants::spin_adapt(const int spin, bitset<nbit__> alpha, bitset<nbit__> beta) const {
  if (spin < 0)
    swap(alpha, beta);

  vector<tuple<int, int, int>> out;

  // bit pattern for doubly occupied orbitals
  bitset<nbit__> common = (alpha & beta);

  bitset<nbit__> alpha_without_common = alpha ^ common;
  bitset<nbit__> beta_without_common = beta ^ common;

  // alpha pattern without highest spin orbitals
  vector<int> salpha_array = bit_to_numbers(alpha_without_common);
  vector<int> ualpha_array;
  if (salpha_array.size() < abs(spin)) throw logic_error("Something is wrong? Determinants::spin_adapt");
  for (int i = 0; i != abs(spin); ++i) {
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

    if (spin >= 0) {
      out.push_back(make_tuple(lexical<1>(ibeta), lexical<0>(ialpha), sign));
    } else {
      out.push_back(make_tuple(lexical<1>(ialpha), lexical<0>(ibeta), sign));
    }
    ++icnt;
  } while (boost::next_combination(open.begin(), open.begin()+nalpha, open.end()));

  // scale to make the vector normalized
  const double factor = 1.0/sqrt(static_cast<double>(icnt));
  return make_pair(out, factor);
}
