//
// BAGEL - Parallel electron correlation program.
// Filename: ras/apply_block.h
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

#ifndef __SRC_RAS_APPLY_BLOCK_H
#define __SRC_RAS_APPLY_BLOCK_H

#include <bitset>
#include <src/util/constants.h>
#include <src/math/algo.h>
#include <src/ciutil/citraits.h>

namespace bagel {

// helper classes for the apply function (this way the main code can be used elsewhere)
namespace RAS {

class Apply_block {
  protected:
    const int orbital_;
    const bool action_;
    const bool spin_;

  private:
    // takes a potential target string, returns whether it is valid and flips the proper bit in the process
    bool test_and_set(std::bitset<nbit__>& bit) {
      bool out = action_ ? bit[orbital_] : !bit[orbital_];
      action_ ? bit.reset(orbital_) : bit.set(orbital_);
      return out;
    }

    int sign(std::bitset<nbit__> bit) const {
      return bagel::sign(bit, orbital_);
    }

  public:
    Apply_block(const int orb, const bool action, const bool spin) : orbital_(orb), action_(action), spin_(spin) {}

    template <class BlockTypeA, class BlockTypeB,
              // checks if A and B are CIBlocks and if they are based on the same StringType
              class = typename std::enable_if<    is_ciblock<BlockTypeA>::value
                                              and is_ciblock<BlockTypeB>::value
                                              and std::is_same<typename BlockTypeA::string_type, typename BlockTypeB::string_type>::value
                                              and std::is_same<typename BlockTypeA::data_type, typename BlockTypeB::data_type>::value
                                             >::type
             >
    void operator()(std::shared_ptr<const BlockTypeA> source, std::shared_ptr<BlockTypeB> target) {
      typedef typename BlockTypeA::data_type DataType;
      if (spin_) {
        const size_t lb = source->lenb();
        assert(lb == target->lenb());
        for (size_t ia = 0; ia < target->lena(); ++ia) {
          std::bitset<nbit__> tabit = target->string_bits_a(ia);
          std::bitset<nbit__> sabit = tabit;
          if (test_and_set(sabit)) {
            const DataType* sourcedata = source->data() + source->stringsa()->lexical_zero(sabit) * lb;
            DataType* targetdata = target->data() + ia * lb;
            const DataType sign = static_cast<DataType>(this->sign(sabit));
            blas::ax_plus_y_n(sign, sourcedata, lb, targetdata);
          }
        }
      }
      else {
        const size_t la = source->lena();
        assert( la == target->lena() );

        const size_t tlb = target->lenb();
        const size_t slb = source->lenb();

        // phase from alpha electrons
        const int alpha_phase = 1 - ((source->stringsa()->nele() & 1) << 1);

        for (auto& tbbit : *target->stringsb()) {
          std::bitset<nbit__> sbbit = tbbit;
          if (test_and_set(sbbit)) {
            DataType* targetdata_base = target->data() + target->stringsb()->lexical_zero(tbbit);
            const DataType* sourcedata_base = source->data() + source->stringsb()->lexical_zero(sbbit);
            const DataType sign = static_cast<DataType>(this->sign(sbbit) * alpha_phase);
            for (size_t i = 0; i < la; ++i) {
              targetdata_base[i*tlb] += sourcedata_base[i*slb] * sign;
            }
          }
        }
      }
    }
};

}
}

#endif
