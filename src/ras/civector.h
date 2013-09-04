//
// BAGEL - Parallel electron correlation program.
// Filename: ras/civector.h
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


#ifndef __BAGEL_RAS_RASCIVECTOR_H
#define __BAGEL_RAS_RASCIVECTOR_H

#include <memory>
#include <tuple>
#include <vector>
#include <list>
#include <map>
#include <bitset>
#include <algorithm>
#include <cassert>

#include <src/ras/determinants.h>
#include <src/util/constants.h>

namespace bagel {

// Contains all the information for a sub block of the CI coefficient matrix
// but does NOT own the data
template <typename DataType>
class RASBlock {
  protected:
    std::shared_ptr<const StringSpace> astrings_;
    std::shared_ptr<const StringSpace> bstrings_;

    DataType* const data_ptr_;

    const int lena_;
    const int lenb_;

    const int offset_;

  public:
    RASBlock(std::shared_ptr<const StringSpace> astrings, std::shared_ptr<const StringSpace> bstrings, DataType* const data_ptr, const int o) :
      astrings_(astrings), bstrings_(bstrings), data_ptr_(data_ptr), lena_(astrings->size()), lenb_(bstrings->size()), offset_(o)
    { }

    const int size() const { return lena_ * lenb_; }
    const int lena() const { return lena_; }
    const int lenb() const { return lenb_; }

    DataType* data() { return data_ptr_; }
    const DataType* data() const { return data_ptr_; }

    DataType& element(const int i) { return data_ptr_[i]; }
    const DataType& element(const int i) const { return data_ptr_[i]; }

    const int index(const std::bitset<nbit__> bbit, const std::bitset<nbit__> abit) const
      { return bstrings_->lexical<0>(bbit) + astrings_->lexical<0>(abit) * lenb_; }
    const int gindex(const std::bitset<nbit__> bbit, const std::bitset<nbit__> abit) const // global index
      { return bstrings_->lexical<0>(bbit) + astrings_->lexical<0>(abit) * lenb_ + offset_; }

    DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) { return element( index(bstring, astring) ); }
    const DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const { return element( index(bstring, astring) ); }

    std::shared_ptr<const StringSpace> stringa() const { return astrings_; }
    std::shared_ptr<const StringSpace> stringb() const { return bstrings_; }
};

template <typename DataType>
class RASCivector {
  using RBlock = RASBlock<DataType>;
  protected:
    std::unique_ptr<DataType[]> data_;
    std::vector<std::shared_ptr<RBlock>> blocks_;

    std::shared_ptr<const RASDeterminants> det_;

    const size_t size_;

    const int hpaddress(const int na, const int nb) const {
      const int N = na + nb;
      return ( (N*(N+1))/2 + nb );
    }

  public:
    RASCivector(std::shared_ptr<const RASDeterminants> det) : det_(det), size_(det->size()) {
      data_ = std::unique_ptr<DataType[]>(new DataType[size_]);
      std::fill_n(data_.get(), size_, 0.0);

      DataType* cc_ptr = data_.get();
      int sz = 0;
      for (auto& ipair : det->stringpairs()) {
        if ( ipair.first && ipair.second ) {
           blocks_.push_back(std::make_shared<RBlock>(ipair.first, ipair.second, cc_ptr, sz));
           cc_ptr += blocks_.back()->size();
           sz += blocks_.back()->size();
        }
        else {
          blocks_.push_back(std::shared_ptr<RBlock>());
        }
      }
    }

    RASCivector(const RASCivector<DataType>& o) : RASCivector(o.det_) {
      std::copy_n(o.data(), size_, data_.get());
    }

    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }

    const int index(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block(bstring, astring)->gindex(bstring, astring);
    }

    DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block(bstring, astring)->element(bstring, astring);
    }
    const DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block(bstring, astring)->element(bstring, astring);
    }

    const std::vector<std::shared_ptr<RBlock>>& blocks() const { return blocks_; }
    std::vector<std::shared_ptr<RBlock>>& blocks() { return blocks_; }

    std::shared_ptr<RBlock> block(const int nha, const int nhb, const int npa, const int npb) {
      const int lp = det_->lenparts();
      return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
    }
    std::shared_ptr<RBlock> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<RBlock> block(std::shared_ptr<const StringSpace> beta, std::shared_ptr<const StringSpace> alpha) {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    std::shared_ptr<const RBlock> block(const int nha, const int nhb, const int npa, const int npb) const {
      const int lp = det_->lenparts();
      return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
    }
    std::shared_ptr<const RBlock> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<const RBlock> block(std::shared_ptr<const StringSpace> beta, std::shared_ptr<const StringSpace> alpha) const {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    template <int spin>
    const std::vector<std::shared_ptr<RBlock>> allowed_blocks(const std::bitset<nbit__> bit) { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }

    template <int spin>
    const std::vector<std::shared_ptr<RBlock>> allowed_blocks(const int nh, const int np) {
      std::vector<std::shared_ptr<RBlock>> out;
      for (int jp = 0; jp + np <= det_->max_particles(); ++jp) {
        for (int ih = 0; ih + nh <= det_->max_holes(); ++ih) {
          std::shared_ptr<RBlock> blk;
          if (spin == 0) blk = block(nh, ih, np, jp);
          else           blk = block(ih, nh, jp, np);

          if (blk) out.push_back(blk);
        }
      }
      return out;
    }

    const size_t size() const { return size_; }
    void zero() { std::fill_n(data_.get(), size_, 0.0); }

    std::shared_ptr<const RASDeterminants> det() const { return det_; }
    std::shared_ptr<RASCivector<DataType>> clone() const { return std::make_shared<RASCivector<DataType>>(det_); }
    std::shared_ptr<RASCivector<DataType>> transpose(std::shared_ptr<const RASDeterminants> det = std::shared_ptr<const RASDeterminants>()) const {
      if (!det) det = det_->transpose();
      auto out = std::make_shared<RASCivector<DataType>>(det);
      for (auto& iblock : blocks_) {
        if (iblock)
          mytranspose_(iblock->data(), iblock->lenb(), iblock->lena(), out->block(iblock->stringa(), iblock->stringb())->data(), 1.0);
      }
      return out;
    }

    DataType dot_product(const RASCivector<DataType>& o) const {
      assert( (*det_) == (*o.det_) );
      return std::inner_product(data_.get(), data_.get() + size_, o.data(), 0.0);
    }

    double norm() const { return std::sqrt(dot_product(*this)); }
    double variance() const { return norm() / size_; }

    void scale(const DataType a) { std::transform( data(), data() + size_, data(), [&a] (DataType p) { return a * p; } ); }
    void ax_plus_y(const DataType a, const RASCivector<DataType>& o)
      { std::transform( o.data(), o.data() + size_, data(), data(), [&a] (DataType p, DataType q) { return (a*p + q); } ); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    double spin_expectation() const { assert(false); return 0.0; } // returns < S^2 >
    std::shared_ptr<RASCivector<DataType>> spin() const { assert(false); return std::shared_ptr<RASCivector<DataType>>();} // returns S^2 | civec >
    std::shared_ptr<RASCivector<DataType>> spin_lower(std::shared_ptr<const RASDeterminants> target_det = std::shared_ptr<RASDeterminants>()) const
      { assert(false); return std::shared_ptr<RASCivector<DataType>>(); } // S_-
    std::shared_ptr<RASCivector<DataType>> spin_raise(std::shared_ptr<const RASDeterminants> target_det = std::shared_ptr<RASDeterminants>()) const
      { assert(false); return std::shared_ptr<RASCivector<DataType>>(); } // S_+
    void spin_decontaminate(const double thresh = 1.0e-8) { assert(false); }

    void project_out(std::shared_ptr<const RASCivector<DataType>> o) { ax_plus_y(-dot_product(*o), *o); }

    double orthog(std::list<std::shared_ptr<const RASCivector<DataType>>> c) {
      for (auto& iter : c)
        project_out(iter);
      const double norm = this->norm();
      const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
      scale(DataType(scal));
      return norm;
    }

    double orthog(std::shared_ptr<const RASCivector<DataType>> o) {
      return orthog(std::list<std::shared_ptr<const RASCivector<DataType>>>{o});
    }

    void print(const double thr) const {
      const DataType* i = data();
      // multimap sorts elements so that they will be shown in the descending order in magnitude
      std::multimap<double, std::tuple<DataType, std::bitset<nbit__>, std::bitset<nbit__>>> tmp;
      for (auto& iblock : blocks_) {
        if (!iblock) continue;
        for (auto& ia : *iblock->stringa()) {
          for (auto& ib : *iblock->stringb()) {
            if (std::abs(*i) > thr)
              tmp.insert(std::make_pair(-std::abs(*i), std::make_tuple(*i, ia, ib)));
            ++i;
          }
        }
      }
      for (auto& iter : tmp)
        std::cout << "       " << det_->print_bit(std::get<1>(iter.second), std::get<2>(iter.second))
                  << "  " << std::setprecision(10) << std::setw(15) << std::get<0>(iter.second) << std::endl;
    }
};

template<> double RASCivector<double>::spin_expectation() const; // returns < S^2 >
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin() const; // returns S^2 | civec >
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(std::shared_ptr<const RASDeterminants>) const; // S_-
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(std::shared_ptr<const RASDeterminants>) const; // S_+
template<> void RASCivector<double>::spin_decontaminate(const double thresh);

using RASCivec = RASCivector<double>;
// RASZCivec may come at some later time

}

#endif
