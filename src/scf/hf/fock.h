//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fock.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __BAGEL_SRC_SCF_FOCK_H
#define __BAGEL_SRC_SCF_FOCK_H

#include <src/df/df.h>
#include <src/integral/libint/libint.h>
#include <src/integral/rys/eribatch.h>
#include <src/scf/hf/fock_base.h>

namespace bagel {

template<int DF>
class Fock : public Fock_base {
  protected:
    void fock_two_electron_part(std::shared_ptr<const Matrix> den = nullptr);
    void fock_two_electron_part_with_coeff(const MatView coeff, const bool rhf, const double scale_ex, const double scale_coulomb);

    // when DF gradients are requested
    bool store_half_;
    std::shared_ptr<DFHalfDist> half_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Fock_base>(*this) & store_half_;
    }

  public:
    Fock() { }
    // Fock operator for DF cases
    template<int DF1 = DF, class = typename std::enable_if<DF1==1>::type>
    Fock(std::shared_ptr<const Geometry> a, std::shared_ptr<const Matrix> prev, std::shared_ptr<const Matrix> den,
         const MatView ocoeff, const bool store = false, const bool rhf = false, const double scale_ex = 1.0, const double scale_coulomb = 1.0)
     : Fock_base(a,prev,den), store_half_(store) {
      fock_two_electron_part_with_coeff(ocoeff, rhf, scale_ex, scale_coulomb);
      fock_one_electron_part();
    }
    // the same as above.
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    Fock(std::shared_ptr<const Geometry> a, std::shared_ptr<const Matrix> prev, std::shared_ptr<const Matrix> den,
         std::shared_ptr<T> ocoeff, const bool store = false, const bool rhf = false, const double scale_ex = 1.0, const double scale_coulomb = 1.0)
     : Fock(a,prev,den,*ocoeff,store,rhf,scale_ex,scale_coulomb) {
    }

    // Fock operator
    template<int DF1 = DF, class = typename std::enable_if<DF1==1 or DF1==0>::type>
    Fock(std::shared_ptr<const Geometry> a, std::shared_ptr<const Matrix> prev, std::shared_ptr<const Matrix> den, const std::vector<double>& d)
     : Fock(a,prev,den,den,d) {
    }

    // Fock operator with a different density matrix for exchange
    template<int DF1 = DF, class = typename std::enable_if<DF1==1 or DF1==0>::type>
    Fock(std::shared_ptr<const Geometry> a, std::shared_ptr<const Matrix> prev, std::shared_ptr<const Matrix> den, std::shared_ptr<const Matrix> ex, const std::vector<double>& d)
     : Fock_base(a,prev,den,d), store_half_(false) {
      fock_two_electron_part(ex);
      fock_one_electron_part();
    }

    std::shared_ptr<DFHalfDist> half() const { return half_; }
};

// specialized for non-DF cases
template <>
void Fock<0>::fock_two_electron_part(std::shared_ptr<const Matrix> den);

}

extern template class bagel::Fock<0>;
extern template class bagel::Fock<1>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Fock<0>)
BOOST_CLASS_EXPORT_KEY(bagel::Fock<1>)

#endif
