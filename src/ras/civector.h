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

//#include <iostream>
//#include <iomanip>
#include <memory>
#include <tuple>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
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

  public:
    RASBlock(std::shared_ptr<const StringSpace> astrings, std::shared_ptr<const StringSpace> bstrings, DataType* const data_ptr) :
      astrings_(astrings), bstrings_(bstrings), data_ptr_(data_ptr), lena_(astrings->size()), lenb_(bstrings->size())
    { }

    const int size() const { return lena_ * lenb_; }
    const int lena() const { return lena_; }
    const int lenb() const { return lenb_; }

    DataType* data() { return data_ptr_; }
    const DataType* data() const { return data_ptr_; }

    DataType& element(const int i) { return data_ptr_[i]; }
    const DataType& element(const int i) const { return data_ptr_[i]; }

    DataType& element(const std::bitset<nbit__> astring, const std::bitset<nbit__> bstring) {
      return element( bstrings_->lexical<0>(bstring) + astrings_->lexical<0>(astring) * lenb_ ); }
    const DataType& element(const std::bitset<nbit__> astring, const std::bitset<nbit__> bstring) const {
      return element( bstrings_->lexical<0>(bstring) + astrings_->lexical<0>(astring) * lenb_ ); }

    std::shared_ptr<const StringSpace> stringa() const { return astrings_; }
    std::shared_ptr<const StringSpace> stringb() const { return bstrings_; }
};

using RBlock = RASBlock<double>;

template <typename DataType>
class RASCivector {
  protected:
    std::unique_ptr<DataType[]> data_;

    std::vector<RASBlock<DataType>> blocks_;

    const size_t size_;

    const int hpaddress(const int na, const int nb) {
      const int N = na + nb;
      return ( (N*(N+1))/2 + nb );
    }

  public:
    RASCivector(std::shared_ptr<const RASDeterminants> det) : size_(det->size()) {
      data_ = std::unique_ptr<DataType[]>(new DataType[size_]);
      std::fill_n(data_.get(), size_, 0.0);

      DataType* cc_ptr = data_.get();
      for (auto& ipair : det->stringpairs()) {
        blocks_.emplace_back(ipair.first, ipair.second, cc_ptr); cc_ptr += blocks_.back().size();
      }
    }

    DataType* data() const { return data_.get(); }
    const DataType* data() const { return data_.get(); }

    DataType& element(const std::bitset<nbit__> astring, const std::bitset<nbit__> bstring) {
      return block(astring, bstring).element(astring, bstring);
    }
    const DataType& element(const std::bitset<nbit__> astring, const std::bitset<nbit__> bstring) const {
      return block(astring, bstring).element(astring, bstring);
    }

    std::vector<RASBlock>& blocks() const { return blocks_; }

    RASBlock<DataType>& block(const int nha, const int nhb, const int npa, const int npb) {
      const int lp = det_->lenparts();
      return block( hpaddress(npa, npb) + lp * hpaddress(nha, nhb) );
    }
    RASBlock<DataType>& block(const std::bitset<nbit__> astring, const std::bitset<nbit__> bstring) {
      return block( det_->nholes(astring), det_->nparticles(astring), det_->nholes(bstring), det_->nparticles(bstring) );
    }
    RASBlock<DataType>& block(std::shared_ptr<const StringSpace> alpha, std::shared_ptr<const StringSpace> beta) {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    template <int spin>
    const std::vector<RASBlock<DataType>> allowed_blocks(const std::bitset<nbit__> bit) { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }

    template <int spin>
    const std::vector<RASBlock<DataType>> allowed_blocks(const int nh, const int np) {
      std::vector<RASBlock<DataType>> out;
      for (int jp = 0; jp + np < det_->max_particles(); ++p) {
        for (int ih = 0; ih + nh < det_->max_holes(); ++ih) {
          if (spin == 0) out.push_back(block(nh, ih, np, jp));
          else           out.push_back(block(ih, nh, jp, np));
        }
      }
      return out;
    }

    const size_t size() const { return size_; }

    void zero() { std::fill_n(data_.get(), size_, 0.0); }

    std::shared_ptr<const RASDeterminants> det() const { return det_; }
    std::shared_ptr<RASCivector<DataType>> clone() const { return std::make_shared<RASCivector<DataType>>(det_); }
    std::shared_ptr<RASCivector<DataType>> transpose(std::shared_ptr<RASDeterminants> det = std::shared_ptr<RASDeterminants>()) const {
      if (!det) det = det_->transpose();
      auto out = make_shared<RASCivector<DataType>>(det);
      for (auto& iblock : blocks_)
        mytranspose_(iblock.data(), iblock.lenb(), iblock.lena(), out->block(iblock.stringa(), iblock.stringb()).data(), 1.0);

      return out;
    }

    DataType dot_product(const RASCivector<DataType>& o) const {
      assert( (*det_) == (*o.det_) );
      std::inner_product(data_.get(), data_.get() + size_, o->data(), 0.0);
    }

    DataType norm() const { return std::sqrt(std::inner_product(data_.get(), data_.get() + size_, data_.get(), 0.0)); }
    DataType variance() const { return norm() / size_; }

    void scale(const DataType a) { std::transform( data(), data() + size_, data(), [&a] (DataType p) { return a * p; } ); }
    void ax_plus_y(const DataType a, const RASCivector<DataType>& o) { std::transform( o.data(), o.data() + size_, data(), [&a] (DataType p, DataType q) { return a*p + q; } ); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    double spin_expectation() const { assert(false); return 0.0; } // returns < S^2 >
    std::shared_ptr<RASCivector<DataType>> spin() const { assert(false); return std::shared_ptr<RASCivector<DataType>>();} // returns S^2 | civec >
    std::shared_ptr<RASCivector<DataType>> spin_lower(std::shared_ptr<const Determinants> target_det = std::shared_ptr<Determinants>()) const
      { assert(false); return std::shared_ptr<RASCivector<DataType>>(); } // S_-
    std::shared_ptr<RASCivector<DataType>> spin_raise(std::shared_ptr<const Determinants> target_det = std::shared_ptr<Determinants>()) const
      { assert(false); return std::shared_ptr<RASCivector<DataType>>(); } // S_+
    void spin_decontaminate(const double thresh = 1.0e-12) { assert(false); }

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
        for (auto& ia : *iblock.stringa()) {
          for (auto& ib : *iblock.stringb()) {
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
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(std::shared_ptr<const Determinants>) const; // S_-
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(std::shared_ptr<const Determinants>) const; // S_+
template<> void RASCivector<double>::spin_decontaminate(const double thresh);

using RASCivec = RASCivector<double>;
// RASZCivec may come at some later time

}

#endif
