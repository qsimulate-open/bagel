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

#include <list>
#include <map>
#include <cassert>
#include <iostream>
#include <iomanip>

#include <src/ras/dvector_base.h>
#include <src/ras/determinants.h>
#include <src/math/algo.h>

namespace bagel {

// Contains all the information for a sub block of the CI coefficient matrix
// but does NOT own the data
template <typename DataType>
class RASBlock {
  protected:
    std::shared_ptr<const StringSpace> astrings_;
    std::shared_ptr<const StringSpace> bstrings_;

    std::unique_ptr<double[]> data_; // can be empty if Block belongs to a Civector
    DataType* const data_ptr_;

    const size_t offset_;

  public:
    RASBlock(std::shared_ptr<const StringSpace> astrings, std::shared_ptr<const StringSpace> bstrings, DataType* const data_ptr, const size_t o) :
      astrings_(astrings), bstrings_(bstrings), data_ptr_(data_ptr), offset_(o) { }
    RASBlock(std::shared_ptr<const StringSpace> astrings, std::shared_ptr<const StringSpace> bstrings) :
      astrings_(astrings), bstrings_(bstrings), data_(new double[size()]), data_ptr_(data_.get()), offset_(0) { std::fill_n(data(), size(), 0.0); }
    // copy constructor
    RASBlock(RASBlock& o) : astrings_(o.astrings_), bstrings_(o.bstrings_), data_(new double[o.size()]), data_ptr_(data_.get()), offset_(0) {
      std::copy_n(o.data(), size(), data());
    }

    const size_t size() const { return lena() * lenb(); }
    const size_t lena() const { return astrings_->size(); }
    const size_t lenb() const { return bstrings_->size(); }

    DataType* data() { return data_ptr_; }
    const DataType* data() const { return data_ptr_; }

    DataType& element(const size_t i) { return data_ptr_[i]; }
    const DataType& element(const size_t i) const { return data_ptr_[i]; }

    const size_t index(const std::bitset<nbit__> bbit, const std::bitset<nbit__> abit) const
      { return bstrings_->lexical<0>(bbit) + astrings_->lexical<0>(abit) * lenb(); }
    const size_t gindex(const std::bitset<nbit__> bbit, const std::bitset<nbit__> abit) const // global index
      { return bstrings_->lexical<0>(bbit) + astrings_->lexical<0>(abit) * lenb() + offset_; }

    DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) { return element( index(bstring, astring) ); }
    const DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const { return element( index(bstring, astring) ); }

    std::shared_ptr<const StringSpace> stringa() const { return astrings_; }
    std::shared_ptr<const StringSpace> stringb() const { return bstrings_; }
};

// helper classes for the apply function (this way the main code can be used elsewhere)
namespace RAS {
template <typename DataType>
class apply_block_base {
  public: virtual void operator()(std::shared_ptr<const RASBlock<DataType>> source, std::shared_ptr<RASBlock<DataType>> target) = 0;
};

template <typename DataType, bool action, bool spin>
class apply_block_impl : public apply_block_base<DataType> {
  protected:
    const int orbital_;

  private:
    bool condition(std::bitset<nbit__>& bit) {
      bool out = action ? !bit[orbital_] : bit[orbital_];
      action ? bit.set(orbital_) : bit.reset(orbital_);
      return out;
    }

    int sign(std::bitset<nbit__> bit) const {
      static_assert(nbit__ <= sizeof(unsigned long long)*8, "verify Determinants::sign (and other functions)");
      bit &= (1ull << orbital_) - 1ull;
      return (1 - (( bit.count() & 1 ) << 1));
    }

  public:
    apply_block_impl(const int orb) : orbital_(orb) {}
    void operator()(std::shared_ptr<const RASBlock<DataType>> source, std::shared_ptr<RASBlock<DataType>> target) override {
      if(spin) {
          const size_t lb = source->lenb();
          assert(lb == target->lenb());
          const DataType* sourcedata = source->data();
          for (auto& abit : *source->stringa()) {
            std::bitset<nbit__> tabit = abit;
            if (condition(tabit)) { // Also sets bit appropriately
              DataType* targetdata = target->data() + target->stringa()->template lexical<0>(tabit) * lb;
              const DataType sign = static_cast<DataType>(this->sign(abit));
              blas::ax_plus_y_n(sign, sourcedata, lb, targetdata);
            }
            sourcedata += lb;
          }
      }
      else {
        const size_t la = source->lena();
        assert( la == target->lena() );
        const DataType* sourcedata_base = source->data();

        const size_t tlb = target->lenb();
        const size_t slb = source->lenb();

        // phase from alpha electrons
        const int alpha_phase = 1 - ((source->stringa()->nele() & 1) << 1);

        for (auto& bbit : *source->stringb()) {
          const DataType* sourcedata = sourcedata_base;
          std::bitset<nbit__> tbbit = bbit;
          if (condition(tbbit)) {
            DataType* targetdata = target->data() + target->stringb()->template lexical<0>(tbbit);
            const DataType sign = static_cast<DataType>(this->sign(bbit) * alpha_phase);
            for (size_t i = 0; i < la; ++i, targetdata+=tlb, sourcedata+=slb) {
              *targetdata += *sourcedata * sign;
            }
          }
          ++sourcedata_base;
        }
      }
    }
};
}

template <typename DataType>
class RASCivector : public std::enable_shared_from_this<RASCivector<DataType>> {
  public: using DetType = RASDeterminants;
  public: using RBlock = RASBlock<DataType>;
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
      size_t sz = 0;
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

    RASCivector(RASCivector<DataType>&& o) : det_(o.det_), size_(o.size_), data_(std::move(o.data_)),
      blocks_(std::move(o.blocks_)) { }

    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }

    // Copy assignment
    RASCivector<DataType>& operator=(const RASCivector<DataType>& o) {
      assert(*o.det_ == *det_);
      std::copy_n(o.data(), size_, data_.get());
      return *this;
    }

    // Move assignment
    RASCivector<DataType>& operator=(RASCivector<DataType>&& o) {
      assert(*o.det_ == *det_);
      data_ = std::move(o.data_);
      blocks_ = std::move(o.blocks_);
      return *this;
    }

    // Element-wise access. Beware: very slow!
    DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block(bstring, astring)->element(bstring, astring);
    }
    const DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block(bstring, astring)->element(bstring, astring);
    }

    // Access to individual blocks
    const std::vector<std::shared_ptr<RBlock>>& blocks() const { return blocks_; }
    std::vector<std::shared_ptr<RBlock>>& blocks() { return blocks_; }

    std::shared_ptr<RBlock> block(const int nha, const int nhb, const int npa, const int npb) {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        const int lp = det_->lenparts();
        return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
      }
      else return std::shared_ptr<RBlock>();
    }
    std::shared_ptr<RBlock> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<RBlock> block(std::shared_ptr<const StringSpace> beta, std::shared_ptr<const StringSpace> alpha) {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    std::shared_ptr<const RBlock> block(const int nha, const int nhb, const int npa, const int npb) const {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        const int lp = det_->lenparts();
        return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
      }
      else return std::shared_ptr<const RBlock>();
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

    template <int spin>
    const std::vector<std::shared_ptr<const RBlock>> allowed_blocks(const std::bitset<nbit__> bit) const { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }
    template <int spin>
    const std::vector<std::shared_ptr<const RBlock>> allowed_blocks(const std::shared_ptr<const StringSpace> space) const { return allowed_blocks<spin>(space->nholes(), space->nparticles()); }

    template <int spin>
    const std::vector<std::shared_ptr<const RBlock>> allowed_blocks(const int nh, const int np) const {
      std::vector<std::shared_ptr<const RBlock>> out;
      for (int jp = 0; jp + np <= det_->max_particles(); ++jp) {
        for (int ih = 0; ih + nh <= det_->max_holes(); ++ih) {
          std::shared_ptr<const RBlock> blk;
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
          blas::transpose(iblock->data(), iblock->lenb(), iblock->lena(), out->block(iblock->stringa(), iblock->stringb())->data(), 1.0);
      }
      return out;
    }

    // Safe for any structure of blocks.
    DataType dot_product(const RASCivector<DataType>& o) const {
      assert( det_->nelea() == o.det()->nelea() && det_->neleb() == o.det()->neleb() && det_->norb() == o.det()->norb() );
      DataType out(0.0);
      for (auto& iblock : this->blocks()) {
        if (iblock) {
          std::shared_ptr<const RBlock> jblock = o.block(iblock->stringb(), iblock->stringa());
          if (jblock) out += blas::dot_product(iblock->data(), iblock->lena()*iblock->lenb(), jblock->data());
        }
      }
      return out;
    }

    double norm() const { return std::sqrt(dot_product(*this)); }
    double variance() const { return dot_product(*this) / size_; }

    void set_det(std::shared_ptr<const RASDeterminants> det) { det_ = det; }
    void scale(const DataType a) { std::for_each( data(), data() + size_, [&a] (DataType& p) { p *= a; } ); }
    void ax_plus_y(const DataType a, const RASCivector<DataType>& o) { blas::ax_plus_y_n(a, o.data(), size_, data()); }
    void ax_plus_y(const DataType a, std::shared_ptr<const RASCivector<DataType>> o) { ax_plus_y(a, *o); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    double spin_expectation() const { assert(false); return 0.0; } // returns < S^2 >
    std::shared_ptr<RASCivector<DataType>> spin() const { assert(false); return std::shared_ptr<RASCivector<DataType>>();} // returns S^2 | civec >
    std::shared_ptr<RASCivector<DataType>> spin_lower(std::shared_ptr<const RASDeterminants> target_det = std::shared_ptr<RASDeterminants>()) const
      { assert(false); return std::shared_ptr<RASCivector<DataType>>(); } // S_-
    std::shared_ptr<RASCivector<DataType>> spin_raise(std::shared_ptr<const RASDeterminants> target_det = std::shared_ptr<RASDeterminants>()) const
      { assert(false); return std::shared_ptr<RASCivector<DataType>>(); } // S_+
    void spin_decontaminate(const double thresh = 1.0e-8) { assert(false); }

    std::shared_ptr<RASCivector<DataType>> apply(const int orbital, const bool action, const bool spin) const {
      // action: true -> create; false -> annihilate
      // spin: true -> alpha; false -> beta
      std::shared_ptr<const RASDeterminants> sdet = this->det();

      const int ras1 = sdet->ras(0);
      const int ras2 = sdet->ras(1);
      const int ras3 = sdet->ras(2);
      const int norb = sdet->norb();

      // 0 -> RASI, 1 -> RASII, 2 -> RASIII
      const int ras_space = ( orbital >= ras1 ) + (orbital >= ras1 + ras2);

      auto to_array = [] (std::shared_ptr<const RASBlock<DataType>> block) {
        auto sa = block->stringa();
        auto sb = block->stringb();
        return std::array<int, 6>({sa->nholes(), sb->nholes(), sa->nele2(), sb->nele2(), sa->nparticles(), sb->nparticles()});
      };

      auto op_on_array = [&ras_space, &action, &spin] ( std::array<int, 6> in ) {
        const int mod = ( action ? +1 : -1 ) * ( ras_space == 0 ? -1 : 1 );
        std::array<int, 6> out = in;
        out[2*ras_space] += (spin ? mod : 0);
        out[2*ras_space+1] += (spin ? 0 : mod);
        return out;
      };

      std::shared_ptr<RAS::apply_block_base<DataType>> apply_block;
      switch ( 2*static_cast<int>(action) + static_cast<int>(spin) ) {
        case 0:
          apply_block = std::make_shared<RAS::apply_block_impl<DataType, false, false>>(orbital); break;
        case 1:
          apply_block = std::make_shared<RAS::apply_block_impl<DataType, false, true>>(orbital);  break;
        case 2:
          apply_block = std::make_shared<RAS::apply_block_impl<DataType, true, false>>(orbital);  break;
        case 3:
          apply_block = std::make_shared<RAS::apply_block_impl<DataType, true, true>>(orbital);   break;
        default:
          assert(false);
      }

      const int mod = action ? +1 : -1;
      const int telea = sdet->nelea() + ( spin ? mod : 0 );
      const int teleb = sdet->neleb() + ( spin ? 0 : mod );
      const int tholes = std::max(sdet->max_holes() - ( (ras_space == 0) ? mod : 0 ), 0);
      const int tparts = std::max(sdet->max_particles() + ( (ras_space == 2) ? mod : 0), 0);

      auto tdet = std::make_shared<const RASDeterminants>(ras1, ras2, ras3, telea, teleb, tholes, tparts, true);
      auto out = std::make_shared<RASCivector<DataType>>(tdet);

      for (auto& soblock : this->blocks()) {
        if (!soblock) continue;
        std::array<int, 6> tar_array = op_on_array(to_array(soblock));
        if ( std::all_of(tar_array.begin(), tar_array.end(), [] (int i) { return i >= 0; }) ) {
          std::shared_ptr<RASBlock<double>> tarblock = out->block(tar_array[0], tar_array[1], tar_array[4], tar_array[5]);
          if (tarblock) (*apply_block)(soblock, tarblock);
        }
      }

      return out;
    }

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

    void print(const double thr = 0.05) const {
      // multimap sorts elements so that they will be shown in the descending order in magnitude
      std::multimap<double, std::tuple<DataType, std::bitset<nbit__>, std::bitset<nbit__>>> tmp;
      for (auto& iblock : blocks_) {
        if (!iblock) continue;
        double* i = iblock->data();
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

using RASDvec = Dvector_base<RASCivec>;

}

#endif
