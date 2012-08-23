//
// Newint - Parallel electron correlation program.
// Filename: space.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_NEWFCI_SPACE_H
#define __SRC_NEWFCI_SPACE_H

#include <memory>
#include <tuple>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <bitset>
#include <cassert>

#include <src/util/constants.h>

#ifdef USE_SSE42_INTRINSICS
#include <nmmintrin.h>
#endif

/************************************************************************************
*     Note: I've been using as many of the member functions of bitset as possible,  *
*        not necessarily because I think it will be faster (I have no idea) but     *
*        because I just want to. I can always change back to faster routines        *
************************************************************************************/

// implements a space that contains all determinants that can be obtained by adding or removing M electrons from a reference
class Space {
  protected:
    // assuming that the number of active orbitals are the same in alpha and beta. 
    const int norb_;

    const int nelea_; // reference number of alpha electrons
    const int neleb_; // reference number of beta electrons

    const int M_; // number of electrons added or removed from a reference

    const bool compress_;

    std::map<int, std::shared_ptr<NewDeterminants> > detmap_; // For now, all access should be through NewDeterminants objects

    template <int spin>
    void form_link_( std::shared_ptr<NewDeterminants> ndet, std::shared_ptr<NewDeterminants> nplusdet ); // links two NewDeterminants objects

    const int key_(const int a, const int b) { return ( (nelea_ - a)*large__ + (neleb_ - b) ); }
    const int key_(std::shared_ptr<NewDeterminants> det) { return key_((nelea_ - det->nelea()),(neleb_ - det->neleb())); }

    const int sign(std::bitset<nbit__> bit, int i) {
      const std::bitset<nbit__> ii( (1 << (i+1)) - 1 );
      bit = bit & ii;
      return (1 - ((bit.count() & 1 ) << 1));
    }

  public:
    Space(std::shared_ptr<NewDeterminants>, const int M);
    Space(const int norb, const int nelea, const int neleb, const int M, const bool compress = true);
    ~Space() {};

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    int nspin() const { return nelea_ - neleb_; }; 

    std::shared_ptr<NewDeterminants> basedet() { return finddet(nelea_, neleb_); }
    // Caution: This function does not check to make sure i,j is valid
    std::shared_ptr<NewDeterminants> finddet( const int i, const int j ) { auto idet = detmap_.find(key_(i,j)); return idet->second; }

};


/************************************************************************************
*  Template function that forms links between two NewDeterminants objects              *
*     - builds the phiup and phidown lists that connect the given Dets              *
*     - assigns the proper detadd/rem links in the Det objects                      *
************************************************************************************/
template <int spin>
void Space::form_link_( std::shared_ptr<NewDeterminants> ndet, std::shared_ptr<NewDeterminants> nplusdet ) {
  if (spin == Alpha) {
    assert((ndet->nelea()+1 == nplusdet->nelea()) && (ndet->neleb() == nplusdet->neleb()));
  }
  else {
    assert((ndet->neleb()+1 == nplusdet->neleb()) && (ndet->nelea() == nplusdet->nelea()));
  }

  std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > > phiup;
  std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > > phidown;

  /* If space is an issue for these functions, I might be able to reduce the amount used... */
  phiup.resize(norb_);
  int upsize = ( (spin==Alpha) ? nplusdet->lena() : nplusdet->lenb() );
  for (auto iter = phiup.begin(); iter != phiup.end(); ++iter) {
    iter->reserve(upsize);
  }

  phidown.resize(norb_);
  int downsize = ( (spin==Alpha) ? ndet->lena() : ndet->lenb() );
  for (auto iter = phidown.begin(); iter != phidown.end(); ++iter) {
    iter->reserve(downsize);
  }

  std::vector<std::bitset<nbit__> > stringplus = ( (spin==Alpha) ? nplusdet->stringa() : nplusdet->stringb() );
  std::vector<std::bitset<nbit__> > string = ( (spin==Alpha) ? ndet->stringa() : ndet->stringb() );

  for (auto iter = string.begin(); iter != string.end(); ++iter) {
    for (unsigned int i = 0; i != norb_; ++i) { 
      if (!(*iter)[i]) { // creation
        const unsigned int source = ndet->lexical<spin>(*iter); 
        std::bitset<nbit__> nbit = *iter; nbit.set(i); // created.
        const unsigned int target = nplusdet->lexical<spin>(nbit);
        phiup[i].push_back(std::make_tuple(target, sign(nbit, i), source));
      }
    }
  }

  for (auto iter = stringplus.begin(); iter != stringplus.end(); ++iter) {
    for (unsigned int i = 0; i!= norb_; ++i) {
      if ((*iter)[i]) { // annihilation
        const int source = nplusdet->lexical<spin>(*iter);
        std::bitset<nbit__> nbit = *iter; nbit.reset(i); //annihilated.
        const unsigned int target = ndet->lexical<spin>(nbit);
        phidown[i].push_back(std::make_tuple(target, sign(nbit, i), source));
      }
    }
  }

  // finally link
  if (spin == Alpha) {
    nplusdet->detremalpha_ = ndet;
    ndet->detaddalpha_ = nplusdet;
  }
  else {
    nplusdet->detrembeta_ = ndet;
    ndet->detaddbeta_ = nplusdet;
  }
};


#endif
