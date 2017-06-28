//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci.h
// Copyright (C) 2011 Toru Shiozaki
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

// Desc :: The implementation closely follows Knowles and Handy 1984 CPL.
//         It is amazing how easy it is to implement FCI !!
//

#ifndef __BAGEL_FCI_KNOWLESHANDY_H
#define __BAGEL_FCI_KNOWLESHANDY_H

#include <src/ci/fci/fci.h>

namespace bagel {

class KnowlesHandy : public FCI {

  protected:
    // denominator
    void const_denom() override;

    // virtual application of Hamiltonian
    std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec> c, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const override;

    // run-time functions
    void sigma_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_3(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_2b (std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, std::shared_ptr<const MOFile> jop) const;
    void sigma_2c1(std::shared_ptr<Civec> sigma, std::shared_ptr<const Dvec> e) const;
    void sigma_2c2(std::shared_ptr<Civec> sigma, std::shared_ptr<const Dvec> e) const;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<FCI>(*this);
      std::shared_ptr<const Matrix> coeff = jop_->coeff();
      ar << coeff;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<FCI>(*this);
      std::shared_ptr<const Matrix> coeff;
      ar >> coeff;
      update(coeff);
    }

  public:
    KnowlesHandy() { }

    KnowlesHandy(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
        const int ncore = -1, const int nocc = -1, const int nstate = -1, const bool store = false);

    KnowlesHandy(std::shared_ptr<const CIWfn> ci, std::shared_ptr<const Reference> r);

    void update(std::shared_ptr<const Matrix>) override;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::KnowlesHandy)

#endif

