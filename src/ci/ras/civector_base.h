//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector_base.h
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

#ifndef __SRC_RAS_CIVECTOR_BASE_H
#define __SRC_RAS_CIVECTOR_BASE_H

#include <src/ci/ras/determinants.h>

namespace bagel {

// Base class contains logic for block structure of RASCivecs
template <class BlockType,
          class = typename std::enable_if<is_ciblock<BlockType>::value>::type
         >
class RASCivector_base {
  protected:
    std::vector<std::shared_ptr<BlockType>> blocks_;
    std::shared_ptr<const RASDeterminants> det_;

    template <class Func>
    void for_each_block(Func func) { for (auto& i: blocks_) if (i) func(i); }

    template <class Func>
    void for_each_block(Func func) const { for (auto& i: blocks_) if (i) func(i); }

    RASCivector_base(std::shared_ptr<const RASDeterminants> d) : det_(d) {}

  public:
    size_t size() const { return det_->size(); }

    std::shared_ptr<const RASDeterminants> det() const { return det_; }
    void set_det(std::shared_ptr<const RASDeterminants> det) { det_ = det; }

    // Access to vectors of blocks
    const std::vector<std::shared_ptr<BlockType>>& blocks() const { return blocks_; }
    std::vector<std::shared_ptr<BlockType>>& blocks() { return blocks_; }

    // Access to individual blocks
    std::shared_ptr<BlockType> block(const int& nha, const int& nhb, const int& npa, const int& npb) {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        return blocks_[det_->block_address(nha, nhb, npa, npb) ];
      }
      else return nullptr;
    }
    std::shared_ptr<BlockType> block(const std::bitset<nbit__>& bstring, const std::bitset<nbit__>& astring) {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<BlockType> block(const std::shared_ptr<const RASString>& beta, const std::shared_ptr<const RASString>& alpha) {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    std::shared_ptr<const BlockType> block(const int& nha, const int& nhb, const int& npa, const int& npb) const {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        return blocks_[det_->block_address(nha, nhb, npa, npb)];
      }
      else return nullptr;
    }
    std::shared_ptr<const BlockType> block(const std::bitset<nbit__>& bstring, const std::bitset<nbit__>& astring) const {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<const BlockType> block(const std::shared_ptr<const RASString>& beta, const std::shared_ptr<const RASString>& alpha) const {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    // Return set of allowed blocks given an input string or block
    template <int spin>
    const std::vector<std::shared_ptr<BlockType>> allowed_blocks(const std::bitset<nbit__>& bit) { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }

    template <int spin>
    const std::vector<std::shared_ptr<BlockType>> allowed_blocks(const int& nh, const int& np) {
      std::vector<std::shared_ptr<BlockType>> out;
      for (int jp = 0; jp + np <= det_->max_particles(); ++jp) {
        for (int ih = 0; ih + nh <= det_->max_holes(); ++ih) {
          std::shared_ptr<BlockType> blk;
          if (spin == 0) blk = block(nh, ih, np, jp);
          else           blk = block(ih, nh, jp, np);

          if (blk) out.push_back(blk);
        }
      }
      return out;
    }

    template <int spin>
    const std::vector<std::shared_ptr<const BlockType>> allowed_blocks(const std::bitset<nbit__>& bit) const { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }
    template <int spin>
    const std::vector<std::shared_ptr<const BlockType>> allowed_blocks(const std::shared_ptr<const RASString>& space) const { return allowed_blocks<spin>(space->nholes(), space->nparticles()); }

    template <int spin>
    const std::vector<std::shared_ptr<const BlockType>> allowed_blocks(const int& nh, const int& np) const {
      std::vector<std::shared_ptr<const BlockType>> out;
      for (int jp = 0; jp + np <= det_->max_particles(); ++jp) {
        for (int ih = 0; ih + nh <= det_->max_holes(); ++ih) {
          std::shared_ptr<const BlockType> blk;
          if (spin == 0) blk = block(nh, ih, np, jp);
          else           blk = block(ih, nh, jp, np);

          if (blk) out.push_back(blk);
        }
      }
      return out;
    }

};

}

#endif
