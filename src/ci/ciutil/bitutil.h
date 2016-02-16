//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: bitutil.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_CIUTIL_BITUTIL_H
#define __SRC_CIUTIL_BITUTIL_H

#include <bitset>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>

#include <src/util/constants.h>

namespace bagel {
  namespace {

    std::string print_bit(std::bitset<nbit__> bit1, std::bitset<nbit__> bit2, const int start, const int fence) {
      std::string out;
      for (int i = start; i != fence; ++i) {
        if (bit1[i] && bit2[i]) { out += "2"; }
        else if (bit1[i]) { out += "a"; }
        else if (bit2[i]) { out += "b"; }
        else { out += "."; }
      }
      return out;
    }
    std::string print_bit(std::bitset<nbit__> bit1, std::bitset<nbit__> bit2, const int max) {
      return print_bit(bit1, bit2, 0, max);
    }

    std::string print_bit(std::bitset<nbit__> bit, const int start, const int fence) {
      std::string out;
      for (int i = start; i != fence; ++i) { out += bit[i] ? "1" : "."; }
      return out;
    }
    std::string print_bit(std::bitset<nbit__> bit, const int max) {
      return print_bit(bit, 0, max);
    }

    std::vector<int> bit_to_numbers(std::bitset<nbit__> bit) {
      std::vector<int> out;
      for (int i = 0; i != nbit__; ++i) if (bit[i]) out.push_back(i);
      return out;
    }

    std::bitset<nbit__> numbers_to_bit(const std::vector<int>& num) {
      std::bitset<nbit__> out(0);
      for (auto& i : num) out.set(i);
      return out;
    }

    int sign(std::bitset<nbit__> bit, int i, int j) {
      // masking irrelevant bits
      int min, max;
      std::tie(min,max) = std::minmax(i,j);
      bit &= ((~std::bitset<nbit__>(0ull)) << min+1);
      bit &= ((~std::bitset<nbit__>(0ull)) >> (nbit__ - max));
      return 1 - ((bit.count() & 1) << 1);
    }

    int sign(std::bitset<nbit__> bit, int i) {
      bit &= ((~std::bitset<nbit__>(0ull)) >> (nbit__ - i));
      return 1 - ((bit.count() & 1 ) << 1);
    }

  }
}

#endif
