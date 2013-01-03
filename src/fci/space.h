//
// BAGEL - Parallel electron correlation program.
// Filename: space.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_FCI_SPACE_H
#define __SRC_FCI_SPACE_H

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
#include <src/fci/determinants.h>

namespace bagel {

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

    std::map<int, std::shared_ptr<Determinants> > detmap_; // For now, all access should be through Determinants objects

    template <int spin>
    void form_link_( std::shared_ptr<Determinants> ndet, std::shared_ptr<Determinants> nplusdet ); // links two Determinants objects

    const int key_(const int a, const int b) { return ( a*large__ + b ); }
    const int key_(std::shared_ptr<Determinants> det) { return key_(det->nelea() - nelea_, det->neleb() - neleb_); }

    const int sign(std::bitset<nbit__> bit, int i) {
      const std::bitset<nbit__> ii( (1 << (i)) - 1 );
      bit = bit & ii;
      return (1 - ((bit.count() & 1 ) << 1));
    }

  public:
    Space(std::shared_ptr<const Determinants>, const int M);
    Space(const int norb, const int nelea, const int neleb, const int M, const bool compress = true);
    ~Space() {};

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    int nspin() const { return nelea_ - neleb_; }; 

    std::shared_ptr<Determinants> basedet() { return finddet(0, 0); };
    // Caution: This function does not check to make sure i,j is valid
    std::shared_ptr<Determinants> finddet( const int i, const int j ) { auto idet = detmap_.find(key_(i,j)); return idet->second; };

  private:
    void common_init();
};


/************************************************************************************
*  Template function that forms links between two Determinants objects              *
*     - builds the phiup and phidown lists that connect the given Dets              *
*     - assigns the proper detadd/rem links in the Det objects                      *
************************************************************************************/
template <int spin>
void Space::form_link_( std::shared_ptr<Determinants> ndet, std::shared_ptr<Determinants> nplusdet ) {
  if (spin == Alpha) {
    assert((ndet->nelea()+1 == nplusdet->nelea()) && (ndet->neleb() == nplusdet->neleb()));
  }
  else {
    assert((ndet->neleb()+1 == nplusdet->neleb()) && (ndet->nelea() == nplusdet->nelea()));
  }

  std::vector<std::vector<DetMap> > phiup;
  std::vector<std::vector<DetMap> > phidown;

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
        phiup[i].push_back(DetMap(target, sign(nbit, i), source));
      }
    }
  }

  for (auto iter = stringplus.begin(); iter != stringplus.end(); ++iter) {
    for (unsigned int i = 0; i!= norb_; ++i) {
      if ((*iter)[i]) { // annihilation
        const unsigned int source = nplusdet->lexical<spin>(*iter);
        std::bitset<nbit__> nbit = *iter; nbit.reset(i); //annihilated.
        const unsigned int target = ndet->lexical<spin>(nbit);
        phidown[i].push_back(DetMap(target, sign(nbit, i), source));
      }
    }
  }


  // finally link
  if (spin == Alpha) {
    nplusdet->detremalpha_ = ndet;
    nplusdet->phidowna_ = phidown;

    ndet->detaddalpha_ = nplusdet;
    ndet->phiupa_ = phiup;
  }
  else {
    nplusdet->detrembeta_ = ndet;
    nplusdet->phidownb_ = phidown;

    ndet->detaddbeta_ = nplusdet;
    ndet->phiupb_ = phiup;
  }
};

}

#endif
