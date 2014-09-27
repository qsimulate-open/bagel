//
// BAGEL - Parallel electron correlation program.
// Filename: fock_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_LONDON_FOCK_LONDON_H
#define __BAGEL_SRC_LONDON_FOCK_LONDON_H

#include <src/df/complexdf.h>
#include <src/integral/comprys/complexeribatch.h>
#include <src/scf/fock_base.h>

namespace bagel {

template<int DF>
class Fock_London : public Fock_base_London {
  protected:
    void fock_two_electron_part(std::shared_ptr<const ZMatrix> den = nullptr);
    void fock_two_electron_part_with_coeff(const ZMatView coeff, const bool rhf, const double scale_ex);

    // when DF gradients are requested
    bool store_half_;
    std::shared_ptr<ComplexDFHalfDist> half_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Fock_base_London>(*this) & store_half_;
    }

  public:
    Fock_London() { }
    // Fock operator for DF cases
    template<int DF1 = DF, class = typename std::enable_if<DF1==1>::type>
    Fock_London(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c,
                const ZMatView ocoeff, const bool store = false, const bool rhf = false, const double scale_ex = 1.0)
     : Fock_base_London(a,b,c), store_half_(store) {
      fock_two_electron_part_with_coeff(ocoeff, rhf, scale_ex);
      fock_one_electron_part();
    }
    // same as above
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    Fock_London(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c,
                std::shared_ptr<T> o, const bool store = false, const bool rhf = false, const double scale_ex = 1.0) : Fock_London(a,b,c,*o,store,rhf,scale_ex) { }

    // Fock operator
    template<int DF1 = DF, class = typename std::enable_if<DF1==1 or DF1==0>::type>
    Fock_London(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c,
                const std::vector<double>& d) : Fock_London(a,b,c,c,d) {}

    // Fock operator with a different density matrix for exchange
    template<int DF1 = DF, class = typename std::enable_if<DF1==1 or DF1==0>::type>
    Fock_London(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c,
                std::shared_ptr<const ZMatrix> ex, const std::vector<double>& d)
     : Fock_base_London(a,b,c,d), store_half_(false) {
      fock_two_electron_part(ex);
      fock_one_electron_part();
    }

    std::shared_ptr<ComplexDFHalfDist> half() const { return half_; }
};

// specialized for non-DF cases
template <>
void Fock_London<0>::fock_two_electron_part(std::shared_ptr<const ZMatrix> den);

}

extern template class bagel::Fock_London<0>;
extern template class bagel::Fock_London<1>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Fock_London<0>)
BOOST_CLASS_EXPORT_KEY(bagel::Fock_London<1>)

#endif
