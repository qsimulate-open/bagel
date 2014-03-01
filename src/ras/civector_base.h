//
// BAGEL - Parallel electron correlation program.
// Filename: ras/civector_base.h
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

#ifndef __SRC_RAS_CIVECTOR_BASE_H
#define __SRC_RAS_CIVECTOR_BASE_H

#include <src/ciutil/citraits.h>
#include <src/ras/determinants.h>

namespace bagel {

// Base class contains logic for block structure of RASCivecs
template <class BlockType,
          class = typename std::enable_if<is_ciblock<BlockType>::value>::type
         >
class RASCivector_base {
  protected:
    std::vector<std::shared_ptr<BlockType>> blocks_;
    std::shared_ptr<const RASDeterminants> det_;

    const int hpaddress(const int na, const int nb) const {
      const int N = na + nb;
      return ( (N*(N+1))/2 + nb );
    }

    template <class Func>
    void for_each_block(Func func) { for (auto& i: blocks_) if (i) func(i); }

    template <class Func>
    void for_each_block(Func func) const { for (auto& i: blocks_) if (i) func(i); }

    RASCivector_base(std::shared_ptr<const RASDeterminants> d) : det_(d) {}

  public:
    std::shared_ptr<const RASDeterminants> det() const { return det_; }
    void set_det(std::shared_ptr<const RASDeterminants> det) { det_ = det; }

    // Access to vectors of blocks
    const std::vector<std::shared_ptr<BlockType>>& blocks() const { return blocks_; }
    std::vector<std::shared_ptr<BlockType>>& blocks() { return blocks_; }

    // Access to individual blocks
    std::shared_ptr<BlockType> block(const int nha, const int nhb, const int npa, const int npb) {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        const int lp = (det_->max_particles()+1) * (det_->max_particles()+2) / 2;
        return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
      }
      else return std::shared_ptr<BlockType>();
    }
    std::shared_ptr<BlockType> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<BlockType> block(std::shared_ptr<const RASString> beta, std::shared_ptr<const RASString> alpha) {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    std::shared_ptr<const BlockType> block(const int nha, const int nhb, const int npa, const int npb) const {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        const int lp = (det_->max_particles()+1) * (det_->max_particles()+2) / 2;
        return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
      }
      else return std::shared_ptr<const BlockType>();
    }
    std::shared_ptr<const BlockType> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<const BlockType> block(std::shared_ptr<const RASString> beta, std::shared_ptr<const RASString> alpha) const {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    // Return set of allowed blocks given an input string or block
    template <int spin>
    const std::vector<std::shared_ptr<BlockType>> allowed_blocks(const std::bitset<nbit__> bit) { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }

    template <int spin>
    const std::vector<std::shared_ptr<BlockType>> allowed_blocks(const int nh, const int np) {
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
    const std::vector<std::shared_ptr<const BlockType>> allowed_blocks(const std::bitset<nbit__> bit) const { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }
    template <int spin>
    const std::vector<std::shared_ptr<const BlockType>> allowed_blocks(const std::shared_ptr<const RASString> space) const { return allowed_blocks<spin>(space->nholes(), space->nparticles()); }

    template <int spin>
    const std::vector<std::shared_ptr<const BlockType>> allowed_blocks(const int nh, const int np) const {
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
