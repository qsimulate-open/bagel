//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: determinants.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <stdexcept>
#include <iomanip>
#include <src/ci/fci/determinants.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <src/util/math/comb.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Determinants)

using namespace std;
using namespace bagel;

const static Comb comb;

using FCIStringSet = CIStringSet<FCIString>;

Determinants::Determinants(const int norb, const int nelea, const int neleb, const bool compress, const bool mute)
 : Determinants(make_shared<FCIString>(nelea, norb), make_shared<FCIString>(neleb, norb), compress, mute) {
}


Determinants::Determinants(shared_ptr<const FCIString> ast, shared_ptr<const FCIString> bst, const bool compress, const bool mute)
 : Determinants(make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{ast}), make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{bst}), compress, mute) {
}


Determinants::Determinants(shared_ptr<const FCIStringSet> ast, shared_ptr<const FCIStringSet> bst, const bool compress, const bool mute) {
  compress_    = compress;
  alphaspaces_ = ast;
  betaspaces_  = bst;
  size_        = lena() * lenb();

  for (auto& a : *alphaspaces_)
    for (auto& b : *betaspaces_)
      blockinfo_.push_back(make_shared<CIBlockInfo<FCIString>>(a, b));

  phia_ = compress_ ? ast->phi() : ast->uncompressed_phi();
  phia_uncompressed_ = ast->uncompressed_phi();
  phib_ = compress_ ? bst->phi() : bst->uncompressed_phi();
  phib_uncompressed_ = bst->uncompressed_phi();

  if (!mute) {
    const int twoS = abs(nspin());
    const int N = nelea() + neleb();
    const size_t out = (twoS + 1) * comb(norb()+1, (N-twoS)/2) * comb(norb()+1, (norb()-((N+twoS)/2)));
    const size_t ncsfs = out / (norb()+1);

    cout << "  Performs exactly the same way as Knowles & Handy 1984 CPL" << endl << endl;
    cout << "  o alpha-beta strings" << endl;
    cout << "      length: " << setw(13) << lena() + lenb() << endl;
    cout << "  o size of the space " << endl;
    cout << "      determinant space:  " << lena() * lenb() << endl;
    cout << "      spin-adapted space: " << ncsfs << endl << endl;
    cout << "  o single displacement lists (alpha)" << endl;
    cout << "      length: " << setw(13) << phia_->size() << endl;
    cout << "  o single displacement lists (beta)" << endl;
    cout << "      length: " << setw(13) << phib_->size() << endl;
  }
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
  return {out, factor};
}


void Determinants::link(shared_ptr<Determinants> odet, shared_ptr<CIStringSpace<FCIStringSet>> spacea,
                                                       shared_ptr<CIStringSpace<FCIStringSet>> spaceb) {
  shared_ptr<Determinants> plusdet;
  shared_ptr<Determinants> det;

  const bool spinb = this->nelea() == odet->nelea();
  const bool spina = this->neleb() == odet->neleb();
  if (!(spina || spinb)) return; // quick return

  const int de = spina ? this->nelea() - odet->nelea() : this->neleb() - odet->neleb();
  if      (de ==  1) tie(det, plusdet) = make_pair(odet, shared_from_this());
  else if (de == -1) tie(det, plusdet) = make_pair(shared_from_this(), odet);
  else return; // quick return

  // finally link
  if (spina) {
    plusdet->set_remalpha(det);
    plusdet->set_phidowna(spacea->phidown(plusdet->stringspacea()));

    det->set_addalpha(plusdet);
    det->set_phiupa(spacea->phiup(det->stringspacea()));
  } else {
    plusdet->set_rembeta(det);
    plusdet->set_phidownb((nelea()&1) ? spaceb->phidown(plusdet->stringspaceb())->get_minus()
                                      : spaceb->phidown(plusdet->stringspaceb()));
    det->set_addbeta(plusdet);
    det->set_phiupb((nelea()&1) ? spaceb->phiup(det->stringspaceb())->get_minus()
                                : spaceb->phiup(det->stringspaceb()));
  }
}
