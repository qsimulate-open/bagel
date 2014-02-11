//
// BAGEL - Parallel electron correlation program.
// Filename: determinants.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_FCI_DETERMINANTS_H
#define __SRC_FCI_DETERMINANTS_H

#include <src/ciutil/determinants_base.h>

namespace bagel {

class Determinants : public Determinants_base<FCIString>,
                     public std::enable_shared_from_this<Determinants> {
  protected:
    /* Links to other determinant spaces accessible by one annihilation or creation operation */
    std::weak_ptr<Determinants> detaddalpha_;
    std::weak_ptr<Determinants> detaddbeta_;
    std::weak_ptr<Determinants> detremalpha_;
    std::weak_ptr<Determinants> detrembeta_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
    }

  public:
    Determinants() { }
    Determinants(const int norb, const int nelea, const int neleb, const bool compress = true, const bool mute = false);
    Determinants(std::shared_ptr<const FCIString> ast, std::shared_ptr<const FCIString> bst, const bool compress = true, const bool mute = false);
    Determinants(std::shared_ptr<const CIStringSet<FCIString>> ast, std::shared_ptr<const CIStringSet<FCIString>> bst, const bool compress = true, const bool mute = false);
    Determinants(std::shared_ptr<const Determinants> o, const bool compress = true, const bool mute = false) :
      Determinants(o->norb(), o->nelea(), o->neleb(), compress, mute) {} // Shortcut to change compression of Det

    std::shared_ptr<Determinants> clone(const int nelea, const int neleb) const {
      return std::make_shared<Determinants>(norb(), nelea, neleb, false, true);
    }

    std::shared_ptr<Determinants> transpose() const { return std::make_shared<Determinants>(norb(), neleb(), nelea(), compress_, true); }

    std::pair<std::vector<std::tuple<int, int, int>>, double> spin_adapt(const int, std::bitset<nbit__>, std::bitset<nbit__>) const;

    std::shared_ptr<const Determinants> addalpha() const { return detaddalpha_.lock();}
    std::shared_ptr<const Determinants> remalpha() const { return detremalpha_.lock();}
    std::shared_ptr<const Determinants> addbeta() const { return detaddbeta_.lock();}
    std::shared_ptr<const Determinants> rembeta() const { return detrembeta_.lock();}

    template<int Spin>
    size_t lexical(const std::bitset<nbit__>& bit) const { return lexical_zero<Spin>(bit); }

    void link(std::shared_ptr<Determinants> odet, std::shared_ptr<CIStringSpace<CIStringSet<FCIString>>>, std::shared_ptr<CIStringSpace<CIStringSet<FCIString>>>);
    template<int spin> void link(std::shared_ptr<Determinants> odet);
};

template<int spin> void Determinants::link(std::shared_ptr<Determinants> odet) {
  std::shared_ptr<Determinants> plusdet;
  std::shared_ptr<Determinants> det;

  const int de = spin == 0 ? this->nelea() - odet->nelea() : this->neleb() - odet->neleb();
  if      (de ==  1) std::tie(det, plusdet) = make_pair(odet, shared_from_this());
  else if (de == -1) std::tie(det, plusdet) = make_pair(shared_from_this(), odet);
  else throw std::logic_error("Determinants::link failed");

  const int fac = (spin == 1 && (nelea() & 1)) ? -1 : 1;
  CIStringSpace<CIStringSet<FCIString>> space{spin==0?alphaspaces_:betaspaces_, spin==0?odet->alphaspaces_:odet->betaspaces_};
  space.build_linkage(fac);

  // finally link
  if (spin == 0) {
    plusdet->detremalpha_ = det;
    plusdet->phidowna_ = space.phidown(plusdet->alphaspaces_);

    det->detaddalpha_ = plusdet;
    det->phiupa_ = space.phiup(det->alphaspaces_);
  } else {
    plusdet->detrembeta_ = det;
    plusdet->phidownb_ = space.phidown(plusdet->betaspaces_);

    det->detaddbeta_ = plusdet;
    det->phiupb_ = space.phiup(det->betaspaces_);
  }
}

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Determinants)

#endif
