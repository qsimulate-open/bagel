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

#include <stdexcept>
#include <src/newfci/determinants.h>
#include <src/util/comb.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>

#include <bitset>

using namespace std;

static Comb comb;

NewDeterminants::NewDeterminants(const int _norb, const int _nelea, const int _neleb, const bool _compress)
  : norb_(_norb), nelea_(_nelea), neleb_(_neleb), compress_(_compress) {

  bool mute = !compress_;
  if (!mute) cout << "  Performs exactly the same way as Knowles & Handy 1984 CPL" << endl << endl;
  if (!mute) cout << "  o lexical mappings" << endl;
  const_lexical_mapping_();
  if (!mute) cout << "  o alpha-beta strings" << endl;
  const_string_lists_();
  if (!mute) cout << "      length: " << setw(13) << stringa_.size() + stringb_.size() << endl;
  if (!mute) cout << "  o single displacement lists (alpha)" << endl;
  const_phis_<0>(stringa_, phia_);
  if (!mute) cout << "      length: " << setw(13) << phia_.size()*phia_.front().size() << endl;
  if (!mute) cout << "  o single displacement lists (beta)" << endl;
  const_phis_<1>(stringb_, phib_);
  if (!mute) cout << "      length: " << setw(13) << phib_.size()*phib_.front().size() << endl << endl;

}


void NewDeterminants::const_string_lists_() {
  vector<int> data(norb_);
  for (int i=0; i!=norb_; ++i)  data[i] = i;

  const int lengtha = comb.c(norb_, nelea_);
  const int lengthb = comb.c(norb_, neleb_);
  stringa_.resize(lengtha);
  stringb_.resize(lengthb);
  fill(stringa_.begin(), stringa_.end(),bitset<nbit__>(0));
  fill(stringb_.begin(), stringb_.end(),bitset<nbit__>(0));

  vector<bitset<nbit__> >::iterator sa = stringa_.begin(); 
  do {
    for (int i=0; i!=nelea_; ++i) sa->set(data[i]);
    ++sa;
  } while (boost::next_combination(data.begin(), data.begin()+nelea_, data.end()));

  sa = stringb_.begin(); 
  do {
    for (int i=0; i!=neleb_; ++i) sa->set(data[i]);
    ++sa;
  } while (boost::next_combination(data.begin(), data.begin()+neleb_, data.end()));

}


void NewDeterminants::const_lexical_mapping_() {
  // combination numbers up to 31 orbitals (util/comb.h)
  zkl_.resize(nelea_ * norb_ + neleb_ * norb_); 
  fill(zkl_.begin(), zkl_.end(), 0u);

  // this part is 1 offset due to the convention of Knowles & Handy's paper.
  // Just a blind copy from the paper without understanding much, but the code below works. 
  for (int k = 1; k < nelea_; ++k) {
    for (int l = k; l <= norb_-nelea_+k; ++l) {
      for (int m = norb_-l+1; m <= norb_-k; ++m) {
        zkl(k-1, l-1, Alpha) += comb.c(m, nelea_-k) - comb.c(m-1, nelea_-k-1); 
      }
    }
  }
  for (int l = nelea_; l <= norb_; ++l) zkl(nelea_-1, l-1, Alpha) = l - nelea_; 

  if (nelea_ == neleb_) {
    copy(zkl_.begin(), zkl_.begin() + nelea_*norb_, zkl_.begin() + nelea_*norb_);
  } else {
    for (int k = 1; k <= neleb_; ++k)
      for (int l = k; l <= norb_-neleb_+k; ++l)
        for (int m = norb_-l+1; m <= norb_-k; ++m)
          zkl(k-1, l-1, Beta) += comb.c(m, neleb_-k) - comb.c(m-1, neleb_-k-1); 
    for (int l = neleb_; l <= norb_; ++l) zkl(neleb_-1, l-1, Beta) = l - neleb_; 
  }
}


void NewDeterminants::print(const double* const civec, const double thr) const {
  const double* i = civec;
  // multimap sorts elements so that they will be shown in the descending order in magnitude
  multimap<double, tuple<double, bitset<nbit__>, bitset<nbit__> > > tmp;
  for (auto ia = stringa_.begin(); ia != stringa_.end(); ++ia) {
    for (auto ib = stringb_.begin(); ib != stringb_.end(); ++ib, ++i) {
      if (abs(*i) > thr) {
        tmp.insert(make_pair(-abs(*i), make_tuple(*i, *ia, *ib))); 
      }
    }
  }
  for (auto iter = tmp.begin(); iter != tmp.end(); ++iter) {
    cout << "       " << print_bit(get<1>(iter->second), get<2>(iter->second))
         << "  " << setprecision(10) << setw(15) << get<0>(iter->second) << endl; 
  }
}


pair<vector<tuple<int, int, int> >, double> NewDeterminants::spin_adapt(const int spin, const bitset<nbit__> alpha, const bitset<nbit__> beta) const {
  vector<tuple<int, int, int> > out;

  int icnt = 0;

  if (spin == 0) {
    bitset<nbit__> common = (alpha & beta);
    const int nalpha = (alpha^common).count();

    vector<int> open = bit_to_numbers(alpha^beta);
    bitset<nbit__> init_alpha = common;
    for (int i=0; i!=nalpha; ++i) init_alpha &= (1<<open[i]);

    // take a linear combination to make a vector singlet coupled.
    do {
      bitset<nbit__> ialpha = common;
      bitset<nbit__> ibeta = common;
      for (int i=0; i!=nalpha; ++i) ialpha ^= (1<<open[i]);
      for (int i=nalpha; i!=open.size(); ++i) ibeta  ^= (1<<open[i]);
      const double sign = pow(-1.0, ((init_alpha^ialpha) >> 2).count());
      out.push_back(make_tuple(lexical<1>(ibeta), lexical<0>(ialpha), sign));
      ++icnt;
    } while (boost::next_combination(open.begin(), open.begin()+nalpha, open.end()));
  } else {
    throw runtime_error("not yet implemented. NewDeterminants::spin_adapt");
  }

  // scale to make the vector normalized
  const double factor = 1.0/sqrt(static_cast<double>(icnt));
  return make_pair(out, factor);
}

